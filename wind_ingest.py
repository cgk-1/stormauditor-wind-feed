#!/usr/bin/env python3
"""
NOAA URMA wind-gust -> Supabase wind-swath ingester for stormauditor.com.

PATCHED (Hazard Engine v2.3): after computing each date's national daily-max
gust grid, this version also samples that grid at EVERY ASOS/AWOS station
(uncapped, below-floor values included) and stores the values via the
hz_station_bg_ingest RPC (src='ANL'). The Hazard Engine's SAWE-2 estimator
uses these as exact objective-analysis backgrounds and as the interpolated
background on days whose grid peaks fall below the 40 mph storage floor.
Zero extra downloads: the grid is already in memory. If the hz_* RPCs are
absent, the sampling step warns and is skipped; swath ingestion is unchanged.

Everything else is identical to the previous version.

Env (GitHub repo secrets): SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: INGEST_DATE, STATES, HOURS_STEP
Deps: pygrib numpy scipy rasterio shapely requests
"""
import os, json, datetime as dt
from zoneinfo import ZoneInfo
import urllib.request
import numpy as np
import pygrib
import requests
from scipy.interpolate import griddata
from rasterio.features import shapes
from rasterio.transform import from_origin
from shapely.geometry import shape, mapping, Point, MultiPolygon, Polygon
from shapely.prepared import prep
from shapely.ops import unary_union

MS2MPH = 2.2369363
BANDS = [40, 58, 74, 96, 111, 130, 157]
POINT_FLOOR = 40
GRID_RES = 0.02

BASE = "https://noaa-urma-pds.s3.amazonaws.com"

# Dominant IANA timezone per state (same convention as the hail feed, so
# wind and hail explorer dates agree with a homeowner's local clock).
STATE_TZ = {
 "Alabama":"America/Chicago","Arizona":"America/Phoenix","Arkansas":"America/Chicago",
 "California":"America/Los_Angeles","Colorado":"America/Denver","Connecticut":"America/New_York",
 "Delaware":"America/New_York","Florida":"America/New_York","Georgia":"America/New_York",
 "Idaho":"America/Boise","Illinois":"America/Chicago","Indiana":"America/Indiana/Indianapolis",
 "Iowa":"America/Chicago","Kansas":"America/Chicago","Kentucky":"America/New_York",
 "Louisiana":"America/Chicago","Maine":"America/New_York","Maryland":"America/New_York",
 "Massachusetts":"America/New_York","Michigan":"America/Detroit","Minnesota":"America/Chicago",
 "Mississippi":"America/Chicago","Missouri":"America/Chicago","Montana":"America/Denver",
 "Nebraska":"America/Chicago","Nevada":"America/Los_Angeles","New Hampshire":"America/New_York",
 "New Jersey":"America/New_York","New Mexico":"America/Denver","New York":"America/New_York",
 "North Carolina":"America/New_York","North Dakota":"America/Chicago","Ohio":"America/New_York",
 "Oklahoma":"America/Chicago","Oregon":"America/Los_Angeles","Pennsylvania":"America/New_York",
 "Rhode Island":"America/New_York","South Carolina":"America/New_York","South Dakota":"America/Chicago",
 "Tennessee":"America/Chicago","Texas":"America/Chicago","Utah":"America/Denver",
 "Vermont":"America/New_York","Virginia":"America/New_York","Washington":"America/Los_Angeles",
 "West Virginia":"America/New_York","Wisconsin":"America/Chicago","Wyoming":"America/Denver",
}
PERMITTED_STATES = set(STATE_TZ)
BOUNDARY_URL = "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json"
_ALL_STATES_CACHE = None
_GEOM = {}


def load_state_geom(name):
    global _ALL_STATES_CACHE
    if name not in _GEOM:
        if _ALL_STATES_CACHE is None:
            gj = json.loads(urllib.request.urlopen(BOUNDARY_URL, timeout=60).read())
            _ALL_STATES_CACHE = {f["properties"]["name"]: f["geometry"] for f in gj["features"]}
        if name not in _ALL_STATES_CACHE:
            raise RuntimeError(f"no boundary found for {name}")
        _GEOM[name] = shape(_ALL_STATES_CACHE[name]).buffer(0)
    return _GEOM[name]


def gust_slice(date_str, hh):
    """Byte-range download ONLY the GUST message for one hour. Returns m/s array."""
    stem = f"urma2p5.{date_str}/urma2p5.t{hh:02d}z.2dvaranl_ndfd.grb2_wexp"
    try:
        idx = urllib.request.urlopen(f"{BASE}/{stem}.idx", timeout=60).read().decode().splitlines()
    except Exception:
        return None
    start = end = None
    for i, line in enumerate(idx):
        f = line.split(":")
        if len(f) > 3 and f[3] == "GUST":
            start = int(f[1])
            end = int(idx[i + 1].split(":")[1]) - 1 if i + 1 < len(idx) else ""
            break
    if start is None:
        return None
    req = urllib.request.Request(f"{BASE}/{stem}", headers={"Range": f"bytes={start}-{end}"})
    data = None
    for attempt in range(4):
        try:
            data = urllib.request.urlopen(req, timeout=180).read()
            if len(data) >= (end - start) if isinstance(end, int) else True:
                break
        except Exception:
            data = None
        import time as _t; _t.sleep(2 * (attempt + 1))
    if not data:
        return None
    open("/tmp/_g.grib2", "wb").write(data)
    g = pygrib.open("/tmp/_g.grib2"); m = g[1]
    v = m.values.astype("float32")
    if "LATS" not in _GEOM:
        la, lo = m.latlons(); _GEOM["LATS"], _GEOM["LONS"] = la, lo
    g.close()
    return v


DMG_MPH = 40  # damaging-wind threshold for duration counting

_SLICE_CACHE = {}

def gust_slice_cached(ds, hh):
    key = (ds, hh)
    if key not in _SLICE_CACHE:
        _SLICE_CACHE[key] = gust_slice(ds, hh)
    return _SLICE_CACHE[key]


def window_hours_utc(tzname, date_str, step=1):
    """UTC (date_str, hour) pairs covering local calendar day D
    (local midnight -> local midnight, DST-aware) in tzname."""
    tz = ZoneInfo(tzname)
    y, m, d = int(date_str[:4]), int(date_str[4:6]), int(date_str[6:])
    start = dt.datetime(y, m, d, tzinfo=tz).astimezone(dt.timezone.utc)
    end = (dt.datetime(y, m, d, tzinfo=tz) + dt.timedelta(days=1)).astimezone(dt.timezone.utc)
    start = start.replace(minute=0, second=0, microsecond=0)
    hours = []
    cur = start
    while cur < end:
        hours.append((cur.strftime("%Y%m%d"), cur.hour))
        cur += dt.timedelta(hours=step)
    return hours


def window_max_mph(hours, step=1):
    """Max gust (mph) + per-cell hours >= DMG_MPH over the given hour list."""
    dmax = None; dur = None
    for ds, hh in hours:
        v = gust_slice_cached(ds, hh)
        if v is None:
            continue
        if dmax is None:
            dmax = v.copy()
            dur = (v * MS2MPH >= DMG_MPH).astype("int16") * step
        else:
            np.fmax(dmax, v, out=dmax)
            dur += (v * MS2MPH >= DMG_MPH).astype("int16") * step
    if dmax is None:
        return None, None, None, None
    return dmax * MS2MPH, dur, _GEOM.get("LATS"), _GEOM.get("LONS")


def daily_max_mph(date_str, step=1):
    """Legacy 06Z-06Z 'convective day' window — still used for the national
    hazard-engine backgrounds so hz_station_bg / hz_bg_coarse semantics are
    unchanged. Explorer state products now use per-state local-day windows."""
    d0 = dt.datetime.strptime(date_str, "%Y%m%d").date()
    d1 = (d0 + dt.timedelta(days=1)).strftime("%Y%m%d")
    hours = [(date_str, hh) for hh in range(6, 24, step)] + \
            [(d1, hh) for hh in range(0, 6, step)]
    return window_max_mph(hours, step)


def build_state(mph, dur, lats, lons, geom):
    minx, miny, maxx, maxy = geom.bounds
    m = (lats >= miny - 0.2) & (lats <= maxy + 0.2) & (lons >= minx - 0.2) & (lons <= maxx + 0.2)
    if not m.any():
        return [], [], 0, 0
    pts = np.column_stack([lons[m], lats[m]]); vals = mph[m]
    dvals = dur[m] if dur is not None else np.zeros(vals.shape, dtype="int16")
    if float(np.nanmax(vals)) < POINT_FLOOR:
        return [], [], 0, 0

    gx = np.arange(minx, maxx, GRID_RES); gy = np.arange(miny, maxy, GRID_RES)
    GX, GY = np.meshgrid(gx, gy)
    grid = np.nan_to_num(griddata(pts, vals, (GX, GY), method="linear"), nan=0.0)
    cls = np.zeros(grid.shape, dtype=np.int16)
    for i, b in enumerate(BANDS, start=1):
        cls[grid >= b] = i
    cls = np.flipud(cls)
    transform = from_origin(minx, maxy, GRID_RES, GRID_RES)
    band_polys = {}
    for g2, val in shapes(cls.astype("int16"), transform=transform):
        val = int(val)
        if val:
            band_polys.setdefault(val, []).append(shape(g2))
    feats = []
    for val, plist in sorted(band_polys.items()):
        merged = unary_union(plist).intersection(geom).simplify(0.01)
        if merged.is_empty:
            continue
        if merged.geom_type == "Polygon":
            merged = MultiPolygon([merged])
        elif merged.geom_type != "MultiPolygon":
            polys = [g for g in merged.geoms if isinstance(g, Polygon)] if hasattr(merged, "geoms") else []
            if not polys:
                continue
            merged = MultiPolygon(polys)
        feats.append({"band": val, "mph_min": BANDS[val - 1], "geom": mapping(merged)})

    pg = prep(geom); points = []; peak_in = 0.0; dur_in = 0
    for (lon, lat), val, dv in zip(pts, vals, dvals):
        if val >= POINT_FLOOR and pg.contains(Point(float(lon), float(lat))):
            points.append({"lon": round(float(lon), 3), "lat": round(float(lat), 3),
                           "v": int(round(float(val)))})
            if val > peak_in: peak_in = float(val)
            if dv > dur_in: dur_in = int(dv)
    if not points:
        return [], [], 0, 0
    return feats, points, int(round(peak_in)), dur_in


def rpc(base, anon, name, payload):
    import time as _t
    last = ""
    for attempt in range(4):
        try:
            r = requests.post(f"{base}/rest/v1/rpc/{name}",
                              headers={"apikey": anon, "Authorization": f"Bearer {anon}",
                                       "Content-Type": "application/json"},
                              data=json.dumps(payload), timeout=120)
            if r.status_code < 300:
                return r
            last = f"{name} {r.status_code}: {r.text[:200]}"
        except Exception as e:
            last = f"{name} exception: {e}"
        _t.sleep(1.5 * (attempt + 1))
    raise RuntimeError(last)


# ---------------- Hazard Engine v2.3: station background sampling ----------
_HZ_STATIONS = None

def _hz_load_stations(base, anon, secret):
    """Station list from the Hazard Engine (fails soft if RPC absent)."""
    global _HZ_STATIONS
    if _HZ_STATIONS is None:
        try:
            r = rpc(base, anon, "hz_stations_fetch", {"p_secret": secret})
            _HZ_STATIONS = r.json() if r.text and r.text != "null" else []
        except Exception as e:
            print(f"  [warn] hz_stations_fetch unavailable ({e}); "
                  f"station backgrounds skipped")
            _HZ_STATIONS = []
    return _HZ_STATIONS


def sample_station_bg(date_iso, mph, lats, lons, base, anon, secret):
    """Sample the (uncapped) national daily-max gust grid at every station and
    upload as SAWE-2 backgrounds. Zero extra downloads."""
    stations = _hz_load_stations(base, anon, secret)
    if not stations:
        return
    la = np.asarray(lats); lo = np.asarray(lons)
    rows = []
    for st in stations:
        try:
            j = int(np.argmin((la - st["lat"])**2 + (lo - st["lon"])**2))
            yy, xx = np.unravel_index(j, la.shape)
            rows.append({"stid": st["stid"],
                         "bg": int(round(float(mph[yy, xx])))})
        except Exception:
            continue
    try:
        for i in range(0, len(rows), 3000):
            rpc(base, anon, "hz_station_bg_ingest",
                {"p_secret": secret, "p_date": date_iso, "p_src": "ANL",
                 "p_rows": rows[i:i+3000], "p_append": i > 0})
        print(f"  {date_iso}  station backgrounds: {len(rows)} sampled (ANL)")
    except Exception as e:
        print(f"  [warn] station bg upload failed: {e}")
    # v2.6: uncapped coarse field samples of the daily-max analysis
    try:
        sub = mph[::10, ::10]; sla = la[::10, ::10]; slo = lo[::10, ::10]
        yy, xx = (sub >= 5).nonzero()
        pts = [{"lon": round(float(slo[a, b]), 2),
                "lat": round(float(sla[a, b]), 2),
                "v": int(round(float(sub[a, b])))}
               for a, b in zip(yy.tolist(), xx.tolist())]
        for i in range(0, len(pts), 4000):
            rpc(base, anon, "hz_bg_coarse_ingest",
                {"p_secret": secret, "p_date": date_iso, "p_src": "ANL",
                 "p_points": pts[i:i+4000]})
        print(f"  {date_iso}  coarse field samples: {len(pts)} (ANL)")
    except Exception as e:
        print(f"  [warn] coarse upload failed: {e}")
# ---------------------------------------------------------------------------


def process_date(date_str, states, base, anon, secret, step):
    date_iso = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
    _SLICE_CACHE.clear()

    # National hazard-engine backgrounds: unchanged legacy 06Z-06Z window.
    mph_bg, dur_bg, lats, lons = daily_max_mph(date_str, step)
    if mph_bg is None:
        print(f"{date_iso}: no URMA data; skipping.")
        _SLICE_CACHE.clear()
        return 0
    sample_station_bg(date_iso, mph_bg, lats, lons, base, anon, secret)
    del mph_bg, dur_bg
    if os.environ.get("BG_ONLY") == "1":
        _SLICE_CACHE.clear()
        return 1   # backgrounds + field samples written; skip swath rebuild

    # Explorer products: per-state LOCAL calendar day (dominant state timezone),
    # so a state's "July 21" is July 21 on a local clock, matching the hail feed.
    groups = {}
    for st in states:
        if st not in PERMITTED_STATES:
            print(f"  [skip] {st} not permitted")
            continue
        groups.setdefault(STATE_TZ[st], []).append(st)

    stored = 0
    for tzname, group_states in sorted(groups.items()):
        hours = window_hours_utc(tzname, date_str, step)
        mph, dur, lats, lons = window_max_mph(hours, step)
        if mph is None:
            print(f"  {date_iso}  [{tzname}] no URMA data for window; skipping group.")
            continue
        for st in group_states:
            try:
                geom = load_state_geom(st)
                feats, points, peak, dur_hrs = build_state(mph, dur, lats, lons, geom)
                if not feats:
                    continue
                rpc(base, anon, "wind_swath_begin",
                    {"p_secret": secret, "p_state": st, "p_date": date_iso,
                     "p_max_mph": peak, "p_dur_hrs": dur_hrs})
                for feat in feats:
                    rpc(base, anon, "wind_swath_add",
                        {"p_secret": secret, "p_state": st, "p_date": date_iso,
                         "p_feature": feat})
                for i in range(0, len(points), 4000):
                    rpc(base, anon, "ingest_wind_points",
                        {"p_secret": secret, "p_state": st, "p_date": date_iso,
                         "p_points": points[i:i+4000], "p_append": i > 0})
                stored += 1
                print(f"  {date_iso}  {st:16s} [{tzname.split('/')[-1]}] peak {peak:.0f} mph, {len(feats)} band(s)")
                import time as _t; _t.sleep(float(os.environ.get("STATE_PAUSE", "0.4")))
            except Exception as e:
                print(f"  [error] {date_iso} {st}: {e}")
        del mph, dur
    _SLICE_CACHE.clear()
    if stored == 0:
        print(f"{date_iso}: no >= {POINT_FLOOR} mph wind on land in selected state(s).")
    return stored


def main():
    raw = os.environ.get("INGEST_DATE") or \
        (dt.datetime.utcnow().date() - dt.timedelta(days=1)).strftime("%Y%m%d")
    dates = []
    for tok in [d.strip() for d in raw.split(",") if d.strip()]:
        if ":" in tok:
            a, b = tok.split(":")
            d0 = dt.datetime.strptime(a, "%Y%m%d").date()
            d1 = dt.datetime.strptime(b, "%Y%m%d").date()
            cur = d0
            while cur <= d1:
                dates.append(cur.strftime("%Y%m%d"))
                cur += dt.timedelta(days=1)
        else:
            dates.append(tok)
    step = int((os.environ.get("HOURS_STEP") or "1").strip() or "1")
    base = os.environ["SUPABASE_URL"].rstrip("/")
    anon = os.environ["SUPABASE_ANON_KEY"]
    secret = os.environ["INGEST_SECRET"]
    states_env = os.environ.get("STATES")
    states = ([s.strip() for s in states_env.split(",")] if states_env else sorted(PERMITTED_STATES))

    print(f"URMA wind ingest: {len(dates)} date(s), {len(states)} state(s), hour step {step}")
    total = 0
    for d in dates:
        total += process_date(d, states, base, anon, secret, step)

    try:
        rpc(base, anon, "purge_old_wind", {"p_secret": secret})
        print("Purged wind data older than 2 years.")
    except Exception as e:
        print(f"[warn] purge failed: {e}")
    print(f"Done. {total} state-day(s) written across {len(dates)} date(s).")


if __name__ == "__main__":
    main()
