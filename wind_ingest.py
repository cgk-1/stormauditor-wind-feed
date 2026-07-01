#!/usr/bin/env python3
"""
NOAA URMA wind-gust -> Supabase wind-swath ingester for stormauditor.com.

Builds a daily "estimated maximum wind gust" product from NOAA's URMA -- the
Unrestricted Mesoscale Analysis, which NOAA designates the *analysis of record*.
URMA is a 2.5 km, hourly, quality-controlled variational analysis that
assimilates surface station wind/gust observations (ASOS/METAR + mesonets) --
i.e. measured station gusts turned into a grid the industry-standard way. It
covers everyday severe/straight-line wind AND tropical systems in one product.

For a given UTC date it: byte-range downloads ONLY the GUST field from each
hourly analysis (~5-6 MB/hour, not the 75 MB full file), takes the 24-hour max,
classifies into severe + Saffir-Simpson wind bands, resamples the Lambert grid to
regular lat/lon per state, polygonizes -> swaths, clips to the state, and pushes
to Supabase via a locked RPC. It also stores raw grid values as points for
precise per-address lookup, and purges anything older than 2 years.

Runs in GitHub Actions (free), never in Supabase/Lovable (which can't read GRIB2).

Env (GitHub repo secrets):
  SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional:
  INGEST_DATE   YYYYMMDD or comma list (default: yesterday UTC)
  STATES        comma list (default: all permitted)
  HOURS_STEP    hour stride for the daily max (default 1; use 3 for fast backfill)

Deps: pygrib numpy scipy rasterio shapely requests
"""
import os, json, datetime as dt
import urllib.request
import numpy as np
import pygrib
import requests
from scipy.interpolate import griddata
from rasterio.features import shapes
from rasterio.transform import from_origin
from shapely.geometry import shape, mapping, Point
from shapely.prepared import prep
from shapely.ops import unary_union

MS2MPH = 2.2369363
# Band lower bounds in mph: 58 = NWS severe-wind criterion (50 kt); 74/96/111/130/157
# = Saffir-Simpson Cat 1-5. 40 = damaging-wind floor (below it we store nothing).
BANDS = [40, 58, 74, 96, 111, 130, 157]
POINT_FLOOR = 40  # store raw grid values (mph) >= this for precise lookups
GRID_RES = 0.02   # target regular grid resolution for swath polygons (deg)

BASE = "https://noaa-urma-pds.s3.amazonaws.com"

PERMITTED_STATES = {
    "Alabama","Arizona","Arkansas","California","Colorado","Connecticut","Delaware",
    "Florida","Georgia","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky",
    "Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi",
    "Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico",
    "New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania",
    "Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont",
    "Virginia","Washington","West Virginia","Wisconsin","Wyoming",
}
BOUNDARY_URL = "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/{slug}.geojson"
_GEOM = {}


def load_state_geom(name):
    if name not in _GEOM:
        url = BOUNDARY_URL.format(slug=name.lower().replace(" ", ""))
        gj = json.loads(urllib.request.urlopen(url, timeout=60).read())
        _GEOM[name] = shape(gj["geometry"] if "geometry" in gj else gj["features"][0]["geometry"]).buffer(0)
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
    open("/tmp/_g.grib2", "wb").write(urllib.request.urlopen(req, timeout=120).read())
    g = pygrib.open("/tmp/_g.grib2"); m = g[1]
    v = m.values.astype("float32")
    if "LATS" not in _GEOM:
        la, lo = m.latlons(); _GEOM["LATS"], _GEOM["LONS"] = la, lo
    g.close()
    return v


def daily_max_mph(date_str, step=1):
    dmax = None
    for hh in range(0, 24, step):
        v = gust_slice(date_str, hh)
        if v is None:
            continue
        dmax = v if dmax is None else np.fmax(dmax, v)
    if dmax is None:
        return None, None, None
    return dmax * MS2MPH, _GEOM["LATS"], _GEOM["LONS"]


def build_state(mph, lats, lons, geom):
    minx, miny, maxx, maxy = geom.bounds
    m = (lats >= miny - 0.2) & (lats <= maxy + 0.2) & (lons >= minx - 0.2) & (lons <= maxx + 0.2)
    if not m.any():
        return [], [], 0.0
    pts = np.column_stack([lons[m], lats[m]]); vals = mph[m]
    peak = float(np.nanmax(vals))
    if peak < POINT_FLOOR:
        return [], [], int(round(peak))

    # resample Lambert -> regular grid for polygon bands
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
        if not merged.is_empty:
            feats.append({"band": val, "mph_min": BANDS[val - 1], "geom": mapping(merged)})

    # raw points >= floor inside the state for precise lookups
    pg = prep(geom); points = []
    for (lon, lat), val in zip(pts, vals):
        if val >= POINT_FLOOR and pg.contains(Point(float(lon), float(lat))):
            points.append({"lon": round(float(lon), 3), "lat": round(float(lat), 3),
                           "v": int(round(float(val)))})
    return feats, points, int(round(peak))


def rpc(base, anon, name, payload):
    r = requests.post(f"{base}/rest/v1/rpc/{name}",
                      headers={"apikey": anon, "Authorization": f"Bearer {anon}",
                               "Content-Type": "application/json"},
                      data=json.dumps(payload), timeout=90)
    if r.status_code >= 300:
        raise RuntimeError(f"{name} {r.status_code}: {r.text[:200]}")
    return r


def process_date(date_str, states, base, anon, secret, step):
    date_iso = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
    mph, lats, lons = daily_max_mph(date_str, step)
    if mph is None:
        print(f"{date_iso}: no URMA data; skipping.")
        return 0
    stored = 0
    for st in states:
        if st not in PERMITTED_STATES:
            print(f"  [skip] {st} not permitted")
            continue
        try:
            geom = load_state_geom(st)
            feats, points, peak = build_state(mph, lats, lons, geom)
            if not feats:
                continue
            rpc(base, anon, "ingest_wind_swath",
                {"p_secret": secret, "p_state": st, "p_date": date_iso,
                 "p_max_mph": peak, "p_features": feats})
            if points:
                rpc(base, anon, "ingest_wind_points",
                    {"p_secret": secret, "p_state": st, "p_date": date_iso, "p_points": points})
            stored += 1
            print(f"  {date_iso}  {st:16s} peak {peak:.0f} mph, {len(feats)} band(s)")
        except Exception as e:
            print(f"  [error] {date_iso} {st}: {e}")
    if stored == 0:
        print(f"{date_iso}: no >= {POINT_FLOOR} mph wind on land in selected state(s).")
    return stored


def main():
    raw = os.environ.get("INGEST_DATE") or \
        (dt.datetime.utcnow().date() - dt.timedelta(days=1)).strftime("%Y%m%d")
    dates = [d.strip() for d in raw.split(",") if d.strip()]
    step = int(os.environ.get("HOURS_STEP", "1"))
    base = os.environ["SUPABASE_URL"].rstrip("/")
    anon = os.environ["SUPABASE_ANON_KEY"]
    secret = os.environ["INGEST_SECRET"]
    states_env = os.environ.get("STATES")
    states = ([s.strip() for s in states_env.split(",")] if states_env else sorted(PERMITTED_STATES))

    print(f"URMA wind ingest: {len(dates)} date(s), {len(states)} state(s), hour step {step}")
    total = 0
    for d in dates:
        total += process_date(d, states, base, anon, secret, step)

    # 2-year retention: purge anything older than today-2y
    try:
        rpc(base, anon, "purge_old_wind", {"p_secret": secret})
        print("Purged wind data older than 2 years.")
    except Exception as e:
        print(f"[warn] purge failed: {e}")
    print(f"Done. {total} state-day(s) written across {len(dates)} date(s).")


if __name__ == "__main__":
    main()
