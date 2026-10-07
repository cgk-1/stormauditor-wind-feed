#!/usr/bin/env python3
"""
NOAA URMA wind-gust -> Supabase wind-swath ingester for stormauditor.com.

PATCHED (Hazard Engine v2.3): after computing each date's national daily-max
gust grid, this version also samples that grid at EVERY ASOS/AWOS station
(uncapped, below-floor values included) and stores the values via the
hz_station_bg_ingest RPC (src='ANL'). The Hazard Engine's SAWE-2 estimator
uses these as exact objective-analysis backgrounds and as the interpolated
background on days whose grid peaks fall below the 40 mph storage floor.
Zero extra downloads: the grid is already in memory.

STAGE 1 HARDENING (Archive Phase 5, 2026-10-07) - the math is unchanged; the
payloads are byte-identical to the previous version for the same input hours:
  * URMA without .idx files (everything before 2024-04-18): the GRIB file is
    walked message by message with small byte-range reads of each header and
    the gust field is selected by its GRIB metadata (discipline 0 / category 2
    / parameter 22) - exactly the bytes the .idx path range-reads (verified on
    overlap dates). Every gust message is validated: GRIB2, 0/2/22, the pinned
    section-3 grid definition (Lambert 2345x1597, md5 09282c59...), valid
    time = requested hour, no bitmap, finite values in 0-120 m/s, lat/lon bounds.
  * A day is never written from a short window (confirmed bug, 2026-10-07: the
    10:10Z run wrote every day before URMA 04Z/05Z/06Z of D+1 existed - they
    land ~10:58Z/11:57Z/12:58Z - so Central states lost their last local hour,
    Mountain 2, Pacific/Arizona 3, the 06Z-06Z backgrounds 2; after DST ends
    Eastern loses one too). Before writing anything, the run computes the
    exact UTC hours the day needs (every state window via window_hours_utc,
    23/24/25 h on DST days, plus the 06Z-06Z background window) and checks
    each one exists upstream. HOURS_POLICY:
      defer   - (scheduled runs) if any needed hour is not published yet the
                WHOLE day is deferred: nothing written, green run, notice in
                the summary. The zero-state self-heal re-runs it later; the
                first run that has every hour for the newest day is the
                16:10Z dispatch (Pacific/Arizona need t06z, ~12:58Z; t07z,
                ~13:58Z, in winter). A needed hour still missing
                URMA_PUBLISH_GRACE_H (12 h) after its valid time is an error.
      strict  - (explicit DATE runs) any missing hour fails the day loudly.
      legacy  - the pre-2026-10-07 behaviour (short windows), kept ONLY to
                reproduce stored days in parity tests. Never the default.
    The self-heal now ignores obs_only rows (wind_obs_rescue, 12:00Z) and also
    checks stored days for SHORT WINDOWS (written before the last needed URMA
    hour was published, by comparing wind_days.updated_at with the S3
    Last-Modified time). Such days are reported, never silently re-written:
    past days are repaired only with owner sign-off (HEAL_SHORT_WINDOWS_FROM=
    YYYY-MM-DD re-runs short days on/after that date).
  * Per-state / station-background / coarse-sample failures are no longer
    swallowed: complete scopes are written, the job ends non-zero.
  * State boundaries are vendored (data/us-states.json, md5-checked).
  * DRY_RUN=1, DATE=YYYY-MM-DD|START..END, completeness summary and the
    FEED_RESULT json line: see feedguard.py.

DAY_CONVENTION=v4 (Archive Phase 5 Stage 3, 2026-10-07; default v3 = unchanged):
  the explorer products (wind_days / wind_polygons / wind_points) use each
  URMA cell's OWN zone's local day (tzwin.py, data/tz/tz_urma.npz) instead of
  its state's zone (plan T1). The 4-5 zone windows are exactly the hours the
  state groups already download (one shared slice cache): no extra downloads.
  Each cell takes the window max / >=40 mph hour count of its own zone's
  window; the state clip, interpolation, banding, rounding and payload layout
  are unchanged. DST days were already exact (23/25 h windows). The 06Z-06Z
  national backgrounds (hz_station_bg ANL, hz_bg_coarse ANL) are NOT changed
  (plan T6, kept for Phase 5). Test flags (DRY_RUN only): V4_ZONES=state.

Env (GitHub repo secrets): SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: DATE / INGEST_DATE, STATES, HOURS_STEP, HOURS_POLICY, BG_ONLY,
          HEAL_SHORT_WINDOWS_FROM, URMA_PUBLISH_GRACE_H, DAY_CONVENTION (v3|v4),
          V4_ZONES / V4_DST (dry-run test flags),
          DRY_RUN, FEED_OUT_DIR, STATE_PAUSE, HZ_STATIONS_FILE (offline station
          list for dry runs without the secret)
Deps: requirements.txt (exact pins)
"""
import os, json, struct, hashlib, tempfile, time, datetime as dt
from zoneinfo import ZoneInfo
import numpy as np
import pygrib
from scipy.interpolate import griddata
from rasterio.features import shapes
from rasterio.transform import from_origin
from shapely.geometry import shape, mapping, Point, MultiPolygon, Polygon
from shapely.prepared import prep
from shapely.ops import unary_union

import feedguard as fg
import tzwin

MS2MPH = 2.2369363
BANDS = [40, 58, 74, 96, 111, 130, 157]
POINT_FLOOR = 40
GRID_RES = 0.02
UTC = dt.timezone.utc

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
# Vendored 2026-10-07, byte-identical to the runtime download every earlier run used
# (raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json).
HERE = os.path.dirname(os.path.abspath(__file__))
BOUNDARY_FILE = os.path.join(HERE, "data", "us-states.json")
BOUNDARY_MD5 = "56968c4d4db9777511c4fd363e684a5c"
_ALL_STATES_CACHE = None
_GEOM = {}

# URMA 2.5 km gust message contract (identical 2021-10 -> 2026-10).
URMA_SEC3_MD5 = "09282c59d78c78302d6e40f6edd446de"   # Lambert 2345x1597, R 6371200, LoV 265, Dx 2539.703 m
URMA_SHAPE = (1597, 2345)
URMA_LAT = (19.228976, 57.088561)
URMA_LON = (-138.373199, -59.042148)
URMA_MAX_MS = fg.env_float("URMA_MAX_MS", 120.0)
URMA_PUBLISH_GRACE_H = fg.env_float("URMA_PUBLISH_GRACE_H", 12.0)


def load_state_geom(name):
    global _ALL_STATES_CACHE
    if name not in _GEOM:
        if _ALL_STATES_CACHE is None:
            with open(BOUNDARY_FILE, "rb") as fh:
                raw = fh.read()
            got = hashlib.md5(raw).hexdigest()
            if got != BOUNDARY_MD5:
                raise fg.ValidationError(f"state boundary file md5 {got} != pinned {BOUNDARY_MD5}")
            gj = json.loads(raw)
            _ALL_STATES_CACHE = {f["properties"]["name"]: f["geometry"] for f in gj["features"]}
        if name not in _ALL_STATES_CACHE:
            raise RuntimeError(f"no boundary found for {name}")
        _GEOM[name] = shape(_ALL_STATES_CACHE[name]).buffer(0)
    return _GEOM[name]


# ---------------------------------------------------------------- GRIB access
def grib2_head(buf):
    """Parse a GRIB2 message head: discipline, edition, total length, section 3
    bytes, product category/number. None when buf is too short to reach the
    end of section 4's first 11 octets."""
    if len(buf) < 16 or buf[:4] != b"GRIB":
        raise fg.ValidationError(f"not a GRIB message (starts {buf[:4]!r})")
    info = {"discipline": buf[6], "edition": buf[7], "total": struct.unpack(">Q", buf[8:16])[0]}
    p = 16
    while p + 5 <= len(buf):
        L = struct.unpack(">I", buf[p:p + 4])[0]
        sn = buf[p + 4]
        if L < 5:
            raise fg.ValidationError(f"corrupt GRIB section length {L}")
        if sn == 3:
            if p + L > len(buf):
                return None
            info["sec3"] = bytes(buf[p:p + L])
        if sn == 4:
            if p + 11 > len(buf):
                return None
            info["cat"], info["num"] = buf[p + 9], buf[p + 10]
            return info
        p += L
    return None


def check_gust_message(data, label):
    info = grib2_head(data)
    if info is None or "sec3" not in info:
        raise fg.ValidationError(f"{label}: GRIB sections 3/4 not found")
    if info["edition"] != 2 or info["discipline"] != 0 or (info["cat"], info["num"]) != (2, 22):
        raise fg.ValidationError(f"{label}: message is ed{info['edition']} {info['discipline']}/"
                                 f"{info['cat']}/{info['num']}, expected GRIB2 0/2/22 (gust)")
    if info["total"] != len(data) or data[-4:] != b"7777":
        raise fg.ValidationError(f"{label}: truncated message ({len(data)} of {info['total']} bytes)")
    s3 = hashlib.md5(info["sec3"]).hexdigest()
    if s3 != URMA_SEC3_MD5:
        raise fg.ValidationError(f"{label}: grid definition changed (section 3 md5 {s3})")


def _gust_bytes_idx(stem, label):
    """The gust message via the .idx byte offsets; None when no .idx exists."""
    raw = fg.http_get(f"{BASE}/{stem}.idx", timeout=60, retries=4, missing_ok=True,
                      what=f"{label} idx")
    if raw is None:
        return None
    idx = raw.decode().splitlines()
    hits = [i for i, line in enumerate(idx) if len(line.split(":")) > 3 and line.split(":")[3] == "GUST"]
    if len(hits) != 1:
        raise fg.ValidationError(f"{label}: {len(hits)} GUST lines in the .idx (expected 1)")
    i = hits[0]
    start = int(idx[i].split(":")[1])
    end = int(idx[i + 1].split(":")[1]) - 1 if i + 1 < len(idx) else None
    return fg.http_get(f"{BASE}/{stem}", byte_range=(start, end), timeout=180, retries=4,
                       what=f"{label} GUST range")


def _gust_bytes_scan(stem, label):
    """No .idx: walk the file's messages with small range reads of each head,
    select the single 0/2/22 message by its GRIB metadata, range-read it."""
    url = f"{BASE}/{stem}"
    size = int(fg.http_head(url, what=label)["Content-Length"])
    off, found, n = 0, [], 0
    while off < size:
        n += 1
        if n > 200:
            raise fg.ValidationError(f"{label}: more than 200 GRIB messages")
        want = 4096
        while True:
            head = fg.http_get(url, byte_range=(off, min(off + want, size) - 1), timeout=60,
                               retries=4, what=f"{label} head@{off}")
            info = grib2_head(head)
            if info is not None or off + want >= size or want >= 1 << 20:
                break
            want *= 4
        if info is None:
            raise fg.ValidationError(f"{label}: cannot parse the message head at byte {off}")
        if info["total"] < 16 or off + info["total"] > size:
            raise fg.ValidationError(f"{label}: message at {off} claims {info['total']} bytes "
                                     f"(file {size})")
        if (info["discipline"], info["cat"], info["num"]) == (0, 2, 22):
            found.append((off, info["total"]))
        off += info["total"]
    if len(found) != 1:
        raise fg.ValidationError(f"{label}: {len(found)} gust (0/2/22) messages (expected 1)")
    off, total = found[0]
    return fg.http_get(url, byte_range=(off, off + total - 1), timeout=180, retries=4,
                       what=f"{label} GUST range")


HOUR_SOURCE = {}   # (date_str, hh) -> "idx" | "scan"


def gust_slice(date_str, hh):
    """Byte-range download ONLY the GUST message for one hour. Returns the m/s
    array (float32). Raises UpstreamMissing when the hour is not published,
    ValidationError / UpstreamError on anything else."""
    stem = f"urma2p5.{date_str}/urma2p5.t{hh:02d}z.2dvaranl_ndfd.grb2_wexp"
    label = f"URMA {date_str} {hh:02d}Z"
    data = _gust_bytes_idx(stem, label)
    via = "idx"
    if data is None:
        data = _gust_bytes_scan(stem, label)
        via = "scan"
    check_gust_message(data, label)
    fd, path = tempfile.mkstemp(suffix=".grib2")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        g = pygrib.open(path)
        try:
            if g.messages != 1:
                raise fg.ValidationError(f"{label}: {g.messages} messages in the gust slice")
            m = g[1]
            got = (m["discipline"], m["parameterCategory"], m["parameterNumber"], m["gridType"],
                   m["Nx"], m["Ny"], int(m["validityDate"]), int(m["validityTime"]))
            want = (0, 2, 22, "lambert", URMA_SHAPE[1], URMA_SHAPE[0], int(date_str), hh * 100)
            if got != want:
                raise fg.ValidationError(f"{label}: GRIB keys {got} != {want}")
            raw_vals = m.values
            if np.ma.isMaskedArray(raw_vals) and np.ma.is_masked(raw_vals):
                raise fg.ValidationError(f"{label}: gust field has masked points")
            v = raw_vals.astype("float32")
            if "LATS" not in _GEOM:
                la, lo = m.latlons()
                if (la.shape != URMA_SHAPE or abs(la.min() - URMA_LAT[0]) > 1e-3
                        or abs(la.max() - URMA_LAT[1]) > 1e-3 or abs(lo.min() - URMA_LON[0]) > 1e-3
                        or abs(lo.max() - URMA_LON[1]) > 1e-3):
                    raise fg.ValidationError(f"{label}: lat/lon bounds {la.min()},{la.max()},"
                                             f"{lo.min()},{lo.max()} not the URMA grid")
                _GEOM["LATS"], _GEOM["LONS"] = la, lo
        finally:
            g.close()
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass
    if v.shape != URMA_SHAPE:
        raise fg.ValidationError(f"{label}: shape {v.shape}")
    if not np.isfinite(v).all():
        raise fg.ValidationError(f"{label}: {int((~np.isfinite(v)).sum())} non-finite values")
    vmin, vmax = float(v.min()), float(v.max())
    if vmin < 0 or vmax > URMA_MAX_MS:
        raise fg.ValidationError(f"{label}: gust range {vmin}..{vmax} m/s outside 0..{URMA_MAX_MS}")
    HOUR_SOURCE[(date_str, hh)] = via
    return v


DMG_MPH = 40  # damaging-wind threshold for duration counting

_SLICE_CACHE = {}


def gust_slice_cached(ds, hh):
    """Array, or None when the hour is not published (404). Any other problem
    raises (and is cached so it is reported once per run)."""
    key = (ds, hh)
    if key not in _SLICE_CACHE:
        try:
            _SLICE_CACHE[key] = gust_slice(ds, hh)
        except fg.UpstreamMissing:
            _SLICE_CACHE[key] = None
        except fg.FeedError as e:
            _SLICE_CACHE[key] = e
    v = _SLICE_CACHE[key]
    if isinstance(v, Exception):
        raise v
    return v


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


def window_max_mph(hours, step=1, missing=None):
    """Max gust (mph) + per-cell hours >= DMG_MPH over the given hour list.
    Hours that are not published are appended to `missing` (the caller applies
    the hours policy); any other fetch/validation problem raises."""
    dmax = None; dur = None
    for ds, hh in hours:
        v = gust_slice_cached(ds, hh)
        if v is None:
            if missing is not None:
                missing.append((ds, hh))
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


def bg_window_hours(date_str, step=1):
    d0 = dt.datetime.strptime(date_str, "%Y%m%d").date()
    d1 = (d0 + dt.timedelta(days=1)).strftime("%Y%m%d")
    return [(date_str, hh) for hh in range(6, 24, step)] + \
           [(d1, hh) for hh in range(0, 6, step)]


def daily_max_mph(date_str, step=1, missing=None):
    """Legacy 06Z-06Z 'convective day' window — still used for the national
    hazard-engine backgrounds so hz_station_bg / hz_bg_coarse semantics are
    unchanged. Explorer state products now use per-state local-day windows."""
    return window_max_mph(bg_window_hours(date_str, step), step, missing)


def hours_ok(run, key, scope, hours, missing, policy, have_data):
    """Apply the hours policy to one window. True = the window may be written."""
    expected = len(hours)
    got = expected - len(missing)
    run.set_received(key, f"hours {scope}", f"{got}/{expected}")
    vias = sorted({HOUR_SOURCE.get(h) for h in hours if h in HOUR_SOURCE} - {None})
    if vias and vias != ["idx"]:
        run.note(key, f"{scope}: URMA read via {'+'.join(vias)} (no .idx for some hours)")
    if not missing:
        return True
    lst = ", ".join(f"{ds[4:6]}-{ds[6:]} {hh:02d}Z" for ds, hh in missing)
    now = dt.datetime.now(UTC)
    recent = all(now - dt.datetime.strptime(ds, "%Y%m%d").replace(tzinfo=UTC)
                 - dt.timedelta(hours=hh) < dt.timedelta(hours=URMA_PUBLISH_GRACE_H)
                 for ds, hh in missing)
    if policy == "legacy" and recent and have_data:
        run.warn(key, f"{scope}: {len(missing)} of {expected} URMA hours not published yet ({lst}) - "
                      f"window ingested WITHOUT them (HOURS_POLICY=legacy; strict would wait)")
        run.day(key).setdefault("partial_windows", {})[scope] = [f"{ds}{hh:02d}" for ds, hh in missing]
        return True
    run.error(key, scope, f"{len(missing)} of {expected} URMA hours missing ({lst}); "
                          f"window NOT written (HOURS_POLICY={policy})",
              details={"missing": [f"{ds}{hh:02d}" for ds, hh in missing]})
    return False


def needed_hours(date_str, states, step=1, bg=True):
    """Exact set of UTC (date, hour) a local day needs: every state window
    (DST-aware) plus, unless bg=False, the 06Z-06Z background window."""
    need = set(bg_window_hours(date_str, step)) if bg else set()
    for tz in {STATE_TZ[s] for s in states}:
        need.update(window_hours_utc(tz, date_str, step))
    return sorted(need)


def hour_published(ds, hh):
    """(published?, Last-Modified or None) for one URMA hour (HEAD, retried)."""
    url = f"{BASE}/urma2p5.{ds}/urma2p5.t{hh:02d}z.2dvaranl_ndfd.grb2_wexp"
    try:
        h = fg.http_head(url, timeout=30, retries=4, what=f"URMA {ds} {hh:02d}Z")
    except fg.UpstreamMissing:
        return False, None
    lm = h.get("Last-Modified")
    import email.utils
    return True, (email.utils.parsedate_to_datetime(lm) if lm else None)


def _recent(ds, hh):
    valid = dt.datetime.strptime(ds, "%Y%m%d").replace(tzinfo=UTC) + dt.timedelta(hours=hh)
    return dt.datetime.now(UTC) - valid < dt.timedelta(hours=URMA_PUBLISH_GRACE_H)


def preflight(run, key, date_str, states, step, policy, bg=True, need=None):
    """True when every needed hour is published. Otherwise defers the day
    (policy defer, all missing hours recent) or fails it - before ANY write.
    need: explicit hour list (v4); default = the v3 set."""
    need = need if need is not None else needed_hours(date_str, states, step, bg)
    missing = [h for h in need if not hour_published(*h)[0]]
    run.set_received(key, "hours_needed", len(need))
    run.set_received(key, "hours_published", len(need) - len(missing))
    if not missing:
        return True
    lst = ", ".join(f"{ds[4:6]}-{ds[6:]} {hh:02d}Z" for ds, hh in missing)
    tags = [f"{ds}{hh:02d}" for ds, hh in missing]
    if policy == "defer" and all(_recent(*h) for h in missing):
        run.defer(key, f"{len(missing)} of {len(need)} needed URMA hours not published yet ({lst}); "
                       f"nothing written - the 16:10Z dispatch (or the next run) ingests the day",
                  pending=tags)
        return False
    run.error(key, "URMA hours", f"{len(missing)} of {len(need)} needed URMA hours missing ({lst}); "
                                 f"day NOT written (HOURS_POLICY={policy})", details={"missing": tags})
    return False


def stored_window_check(run, key, date_str, states, step):
    """Hour-aware completeness of a STORED day: True when the grid rows were
    written before the last needed URMA hour was published (a short window).
    Reads wind_days.updated_at (feed rows only) + S3 Last-Modified."""
    import requests
    r = requests.get(f"{run.base}/rest/v1/wind_days",
                     params={"valid_date": f"eq.{key}", "obs_only": "is.false",
                             "select": "updated_at", "order": "updated_at.asc", "limit": "1"},
                     headers={"apikey": run.anon, "Authorization": f"Bearer {run.anon}"},
                     timeout=30)
    if r.status_code >= 300:
        raise fg.FeedError(f"wind_days read HTTP {r.status_code}")
    rows = r.json()
    if not rows:
        return False
    written = dt.datetime.fromisoformat(rows[0]["updated_at"].replace("Z", "+00:00"))
    late = []
    for ds, hh in needed_hours(date_str, states, step):
        ok, lm = hour_published(ds, hh)
        if not ok or (lm and lm > written):
            late.append(f"{ds}{hh:02d}")
    if late:
        run.set_received(key, "stored_short_window_hours", late)
    return bool(late)


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


def validate_state_output(state, geom, feats, points, peak, dur_hrs):
    minx, miny, maxx, maxy = geom.bounds
    seen = set()
    for p in points:
        lon, lat, v = p["lon"], p["lat"], p["v"]
        if not (POINT_FLOOR <= v <= 300):
            raise fg.ValidationError(f"{state}: point value {v} mph outside 40-300")
        if not (minx - 0.01 <= lon <= maxx + 0.01 and miny - 0.01 <= lat <= maxy + 0.01):
            raise fg.ValidationError(f"{state}: point {lon},{lat} outside the state bounds")
        if (lon, lat) in seen:
            raise fg.ValidationError(f"{state}: duplicate point {lon},{lat}")
        seen.add((lon, lat))
    if not (POINT_FLOOR <= peak <= 300) or not (0 <= dur_hrs <= 25):
        raise fg.ValidationError(f"{state}: peak {peak} mph / duration {dur_hrs} h out of range")
    last = 0
    for f in feats:
        if not (last < f["band"] <= len(BANDS)) or f["mph_min"] != BANDS[f["band"] - 1]:
            raise fg.ValidationError(f"{state}: unexpected band {f['band']}/{f['mph_min']}")
        if f["geom"]["type"] != "MultiPolygon" or not f["geom"]["coordinates"]:
            raise fg.ValidationError(f"{state}: band {f['band']} geometry {f['geom']['type']}")
        last = f["band"]


# ---------------- Hazard Engine v2.3: station background sampling ----------
_HZ_STATIONS = None


def _hz_load_stations(run):
    """Station list from the Hazard Engine (read-only RPC). A failure is an
    error now (it used to skip the backgrounds silently)."""
    global _HZ_STATIONS
    if _HZ_STATIONS is None:
        path = os.environ.get("HZ_STATIONS_FILE")
        if path:
            with open(path) as fh:
                _HZ_STATIONS = json.load(fh)
        else:
            r = run.rpc_read("hz_stations_fetch", {"p_secret": run.secret})
            _HZ_STATIONS = r.json() if r.text and r.text != "null" else []
        if not isinstance(_HZ_STATIONS, list) or len(_HZ_STATIONS) < 1000:
            n = len(_HZ_STATIONS) if isinstance(_HZ_STATIONS, list) else "?"
            _HZ_STATIONS = None
            raise fg.ValidationError(f"hz_stations_fetch returned {n} stations (expected > 1000)")
    return _HZ_STATIONS


def sample_station_bg(run, key, date_iso, mph, lats, lons):
    """Sample the (uncapped) national daily-max gust grid at every station and
    upload as SAWE-2 backgrounds. Zero extra downloads."""
    stations = _hz_load_stations(run)
    run.expect(key, "stations", len(stations))
    la = np.asarray(lats); lo = np.asarray(lons)
    rows = []
    for st in stations:
        j = int(np.argmin((la - st["lat"])**2 + (lo - st["lon"])**2))
        yy, xx = np.unravel_index(j, la.shape)
        rows.append({"stid": st["stid"],
                     "bg": int(round(float(mph[yy, xx])))})
    bad = [r for r in rows if not (0 <= r["bg"] <= 300)]
    if bad:
        raise fg.ValidationError(f"{len(bad)} station backgrounds out of 0-300 mph, e.g. {bad[:3]}")
    run.set_received(key, "stations", len(rows))
    calls = [("hz_station_bg_ingest",
              {"p_secret": run.secret, "p_date": date_iso, "p_src": "ANL",
               "p_rows": rows[i:i+3000], "p_append": i > 0})
             for i in range(0, len(rows), 3000)]
    run.write(key, "ANL station_bg", calls)
    run.written(key, "ANL station_bg")
    print(f"  {date_iso}  station backgrounds: {len(rows)} sampled (ANL)")
    # v2.6: uncapped coarse field samples of the daily-max analysis
    sub = mph[::10, ::10]; sla = la[::10, ::10]; slo = lo[::10, ::10]
    yy, xx = (sub >= 5).nonzero()
    pts = [{"lon": round(float(slo[a, b]), 2),
            "lat": round(float(sla[a, b]), 2),
            "v": int(round(float(sub[a, b])))}
           for a, b in zip(yy.tolist(), xx.tolist())]
    if not pts or any(not (5 <= p["v"] <= 300) for p in pts):
        raise fg.ValidationError(f"coarse samples: {len(pts)} points or values out of 5-300 mph")
    calls = [("hz_bg_coarse_ingest",
              {"p_secret": run.secret, "p_date": date_iso, "p_src": "ANL",
               "p_points": pts[i:i+4000]})
             for i in range(0, len(pts), 4000)]
    run.write(key, "ANL coarse", calls)
    run.written(key, "ANL coarse")
    print(f"  {date_iso}  coarse field samples: {len(pts)} (ANL)")
# ---------------------------------------------------------------------------


def process_date(run, date_str, states, step, policy, bg_only=False):
    date_iso = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
    key = date_iso
    _SLICE_CACHE.clear()
    run.expect(key, "states", len(states))
    if policy in ("defer", "strict") and not preflight(run, key, date_str, states, step, policy,
                                                       bg=True):
        return 0

    # National hazard-engine backgrounds: unchanged legacy 06Z-06Z window.
    bg_hours = bg_window_hours(date_str, step)
    missing = []
    try:
        mph_bg, dur_bg, lats, lons = daily_max_mph(date_str, step, missing)
        if hours_ok(run, key, "ANL background 06Z-06Z", bg_hours, missing, policy, mph_bg is not None):
            if mph_bg is None:
                raise fg.UpstreamMissing("no URMA hour available for the 06Z-06Z window")
            sample_station_bg(run, key, date_iso, mph_bg, lats, lons)
        del mph_bg, dur_bg
    except Exception as e:
        run.error(key, "ANL background 06Z-06Z", f"{type(e).__name__}: {e}")
    if bg_only:
        _SLICE_CACHE.clear()
        return 0

    # Explorer products: per-state LOCAL calendar day (dominant state timezone),
    # so a state's "July 21" is July 21 on a local clock, matching the hail feed.
    groups = {}
    for st in states:
        groups.setdefault(STATE_TZ[st], []).append(st)
    run.expect(key, "windows", len(groups))

    stored = 0
    for tzname, group_states in sorted(groups.items()):
        scope = f"window {tzname}"
        hours = window_hours_utc(tzname, date_str, step)
        missing = []
        try:
            mph, dur, lats, lons = window_max_mph(hours, step, missing)
            ok = hours_ok(run, key, scope, hours, missing, policy, mph is not None)
            if ok and mph is None:
                raise fg.UpstreamMissing("no URMA hour available for the window")
        except Exception as e:
            run.error(key, scope, f"{type(e).__name__}: {e}")
            ok = False
        if not ok:
            for st in group_states:
                run.day(key)["failed"].append(st)
            continue
        run.receive(key, "windows")
        for st in group_states:
            try:
                geom = load_state_geom(st)
                feats, points, peak, dur_hrs = build_state(mph, dur, lats, lons, geom)
                if not feats:
                    run.empty(key, st)
                    continue
                validate_state_output(st, geom, feats, points, peak, dur_hrs)
                calls = [(SWATH_BEGIN_RPC,
                          {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                           "p_max_mph": peak, "p_dur_hrs": dur_hrs})]
                for feat in feats:
                    calls.append(("wind_swath_add",
                                  {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                                   "p_feature": feat}))
                for i in range(0, len(points), 4000):
                    calls.append(("ingest_wind_points",
                                  {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                                   "p_points": points[i:i+4000], "p_append": i > 0}))
                run.write(key, st, calls)
                run.written(key, st)
                stored += 1
                print(f"  {date_iso}  {st:16s} [{tzname.split('/')[-1]}] peak {peak:.0f} mph, {len(feats)} band(s)")
                if not run.dry_run:
                    time.sleep(float(os.environ.get("STATE_PAUSE", "0.4")))
            except Exception as e:
                run.error(key, st, f"{type(e).__name__}: {e}")
        del mph, dur
    _SLICE_CACHE.clear()
    d = run.day(key)
    if not run.dry_run and d["written"]:
        check_obs_only(run, key, d["written"])
    run.set_received(key, "states_ok", len(d["written"]) + len([s for s in d["empty"] if s in PERMITTED_STATES]))
    if stored == 0 and not d["failed"]:
        print(f"{date_iso}: no >= {POINT_FLOOR} mph wind on land in selected state(s).")
    return stored


# ===================================================================== v4
# DAY_CONVENTION=v4 (Archive Phase 5 Stage 3): own-zone local day per URMA cell
# for the explorer products. Used only when DAY_CONVENTION=v4.

LOWER48_ZONE_IDS = range(1, 22)     # tz MANIFEST ids 1..21 = the lower-48 zones


def needed_hours_v4(date_str, states, step=1, bg=True):
    """v3's hour set plus every lower-48 zone window (identical on CONUS: the
    zone families are the same as the state zones')."""
    need = set(needed_hours(date_str, states, step, bg))
    names = tzwin.zone_names()
    for z in LOWER48_ZONE_IDS:
        need.update(window_hours_utc(names[z], date_str, step))
    return sorted(need)


def process_date_v4(run, date_str, states, step, policy, flags, bg_only=False):
    real = flags["zones"] == "real"
    date_iso = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]}"
    key = date_iso
    _SLICE_CACHE.clear()
    run.expect(key, "states", len(states))
    if policy in ("defer", "strict") and not preflight(
            run, key, date_str, states, step, policy, bg=True,
            need=needed_hours_v4(date_str, states, step, True)):
        return 0

    # National hazard-engine backgrounds: unchanged legacy 06Z-06Z window (T6).
    bg_hours = bg_window_hours(date_str, step)
    missing = []
    try:
        mph_bg, dur_bg, lats, lons = daily_max_mph(date_str, step, missing)
        if hours_ok(run, key, "ANL background 06Z-06Z", bg_hours, missing, policy, mph_bg is not None):
            if mph_bg is None:
                raise fg.UpstreamMissing("no URMA hour available for the 06Z-06Z window")
            sample_station_bg(run, key, date_iso, mph_bg, lats, lons)
        del mph_bg, dur_bg
    except Exception as e:
        run.error(key, "ANL background 06Z-06Z", f"{type(e).__name__}: {e}")
    if bg_only:
        _SLICE_CACHE.clear()
        return 0

    zm = tzwin.zone_map("urma")
    if "LATS" not in _GEOM:      # lat/lon of the grid come with the first decoded slice
        for h in window_hours_utc(STATE_TZ[states[0]], date_str, step):
            try:
                if gust_slice_cached(*h) is not None:
                    break
            except fg.FeedError:
                continue
    if "LATS" not in _GEOM:
        run.error(key, "URMA grid", "no URMA hour could be decoded")
        return 0
    LA, LO = _GEOM["LATS"], _GEOM["LONS"]
    if LA.shape != zm.zone.shape:
        raise fg.ValidationError(f"URMA grid {LA.shape} != zone map {zm.zone.shape}")
    geoms = {st: load_state_geom(st) for st in states}
    st_zone = {st: tzwin.zone_id(STATE_TZ[st]) for st in states}
    st_zones, bbox = {}, {}
    for st in states:
        minx, miny, maxx, maxy = geoms[st].bounds
        m = (LA >= miny - 0.2) & (LA <= maxy + 0.2) & (LO >= minx - 0.2) & (LO <= maxx + 0.2)
        bbox[st] = m
        zs = set(np.unique(zm.zone[m]).tolist()) - {0} if real else set()
        st_zones[st] = zs | {st_zone[st]}
    dg = tzwin.DayGroups(date_str, set().union(*st_zones.values()))
    st_group = {st: tzwin.group_of_window(dg, tzwin.local_window(STATE_TZ[st], date_str)) for st in states}
    st_groups = {st: sorted({int(dg.lut[z]) for z in st_zones[st]}) for st in states}
    used = sorted({g for st in states for g in st_groups[st]})
    run.day(key)["v4_groups"] = dg.describe()
    run.expect(key, "windows", len(used))

    fields, failed_g = {}, {}
    for g in used:
        rep = dg.members[g][0]
        scope = f"window {rep}"
        hours = window_hours_utc(rep, date_str, step)
        missing = []
        try:
            mph, dur, _la, _lo = window_max_mph(hours, step, missing)
            ok = hours_ok(run, key, scope, hours, missing, policy, mph is not None)
            if ok and mph is None:
                raise fg.UpstreamMissing("no URMA hour available for the window")
        except Exception as e:
            run.error(key, scope, f"{type(e).__name__}: {e}")
            ok = False
        if not ok:
            failed_g[g] = scope
            continue
        run.receive(key, "windows")
        fields[g] = (mph, dur)

    if real and fields:
        gmap = dg.lut[zm.zone]
        mph_any = np.max(np.stack([f[0] for f in fields.values()]), axis=0)
        dur_any = np.max(np.stack([f[1] for f in fields.values()]), axis=0)
        # zone-0 cells (>1 deg from US territory) are never stored; they get the
        # max over all windows so the interpolation input is fully defined.
        comp_mph = tzwin.compose(gmap, {g: f[0] for g, f in fields.items()}, mph_any)
        comp_dur = tzwin.compose(gmap, {g: f[1] for g, f in fields.items()}, dur_any)
        del mph_any, dur_any

    # Same iteration (and payload) order as v3: by state tz name, then input order.
    by_tz = {}
    for st in states:
        by_tz.setdefault(STATE_TZ[st], []).append(st)
    stored = 0
    for tzname, group_states in sorted(by_tz.items()):
        for st in group_states:
            bad = [g for g in st_groups[st] if g in failed_g]
            if bad:
                run.day(key)["failed"].append(st)
                continue
            try:
                geom = geoms[st]
                if real:
                    mph, dur = comp_mph, comp_dur
                else:
                    mph, dur = fields[st_group[st]]
                feats, points, peak, dur_hrs = build_state(mph, dur, LA, LO, geom)
                if not feats:
                    run.empty(key, st)
                    continue
                if real:
                    m = bbox[st]
                    zone_at = dict(zip(((round(float(x), 3), round(float(y), 3))
                                        for x, y in zip(LO[m].tolist(), LA[m].tolist())),
                                       zm.zone[m].tolist()))
                    for p in points:
                        z = zone_at.get((p["lon"], p["lat"]))
                        if not z:
                            raise fg.ValidationError(f"{st}: stored cell {p['lon']},{p['lat']} has no US "
                                                     f"zone (zone id {z})")
                validate_state_output(st, geom, feats, points, peak, dur_hrs)
                calls = [(SWATH_BEGIN_RPC,
                          {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                           "p_max_mph": peak, "p_dur_hrs": dur_hrs})]
                for feat in feats:
                    calls.append(("wind_swath_add",
                                  {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                                   "p_feature": feat}))
                for i in range(0, len(points), 4000):
                    calls.append(("ingest_wind_points",
                                  {"p_secret": run.secret, "p_state": st, "p_date": date_iso,
                                   "p_points": points[i:i+4000], "p_append": i > 0}))
                run.write(key, st, calls)
                run.written(key, st)
                stored += 1
                print(f"  {date_iso}  {st:16s} [{tzname.split('/')[-1]}{' +' + str(len(st_groups[st]) - 1) + ' zone window(s)' if len(st_groups[st]) > 1 else ''}] "
                      f"peak {peak:.0f} mph, {len(feats)} band(s)")
                if not run.dry_run:
                    time.sleep(float(os.environ.get("STATE_PAUSE", "0.4")))
            except Exception as e:
                run.error(key, st, f"{type(e).__name__}: {e}")
    _SLICE_CACHE.clear()
    d = run.day(key)
    if not run.dry_run and d["written"]:
        check_obs_only(run, key, d["written"])
    run.set_received(key, "states_ok", len(d["written"]) + len([s for s in d["empty"] if s in PERMITTED_STATES]))
    if stored == 0 and not d["failed"]:
        print(f"{date_iso}: no >= {POINT_FLOOR} mph wind on land in selected state(s).")
    return stored


# Owner-approved 2026-10-07: wind_swath_begin_v2 also clears obs_only when the
# grid lands (a deferred day can be rescued at 12:00Z before its grid exists).
# legacy keeps v1 so test replays reproduce the stored night-one calls.
SWATH_BEGIN_RPC = "wind_swath_begin_v2"


def check_obs_only(run, key, written):
    """wind_obs_rescue (pg_cron 12:00Z) inserts obs_only=true wind_days rows
    for states with no row yet; wind_swath_begin's upsert never resets that
    flag, and the site then hides the grid swath for the day. A deferred day
    can be rescued before its grid lands, so say so loudly."""
    import requests
    try:
        r = requests.get(f"{run.base}/rest/v1/wind_days",
                         params={"valid_date": f"eq.{key}", "obs_only": "is.true", "select": "state"},
                         headers={"apikey": run.anon, "Authorization": f"Bearer {run.anon}"},
                         timeout=30)
        hit = sorted({x["state"] for x in r.json()} & set(written)) if r.status_code < 300 else None
    except Exception as e:
        hit = None
        print(f"  [warn] obs_only check failed: {e}")
    if hit:
        run.warn(key, f"grid written for {hit} but their wind_days rows are still obs_only=true "
                      f"(rescued at 12:00Z before the grid landed): the site hides these swaths "
                      f"until obs_only is reset (needs the DB fix in the Stage 1 report)")


def parse_states():
    states_env = (os.environ.get("STATES") or "").strip()
    if not states_env:
        return sorted(PERMITTED_STATES)
    states = [s.strip() for s in states_env.split(",") if s.strip()]
    bad = [s for s in states if s not in PERMITTED_STATES]
    if bad:
        raise fg.ValidationError(f"STATES has unknown/not permitted state(s): {bad}")
    return states


def main(run):
    step = int((os.environ.get("HOURS_STEP") or "1").strip() or "1")
    bg_only = os.environ.get("BG_ONLY") == "1"
    explicit = fg.requested_dates()
    policy = (os.environ.get("HOURS_POLICY") or "").strip().lower() or \
        ("strict" if explicit is not None else "defer")
    if policy not in ("strict", "defer", "legacy"):
        raise fg.ValidationError(f"HOURS_POLICY={policy!r} (use defer, strict or legacy)")
    global SWATH_BEGIN_RPC
    SWATH_BEGIN_RPC = "wind_swath_begin" if policy == "legacy" else "wind_swath_begin_v2"
    heal_from = (os.environ.get("HEAL_SHORT_WINDOWS_FROM") or "").strip()
    heal_from = fg._one_date(heal_from) if heal_from else None
    run.meta.update({"boundary_md5": BOUNDARY_MD5, "hours_policy": policy, "hours_step": step,
                     "numpy": np.__version__, "pygrib": pygrib.__version__})
    conv = tzwin.convention()
    flags = tzwin.test_flags(run.dry_run) if conv == "v4" else None
    run.meta["day_convention"] = conv
    if conv == "v4":
        run.meta.update({"tzwin_md5": tzwin.module_md5(), "zone_map_md5": tzwin.zone_map("urma").md5,
                         "v4_flags": flags})

    dates = []
    if explicit is not None:
        dates = [d.strftime("%Y%m%d") for d in explicit]
    else:
        # Scheduled run: yesterday, PLUS SELF-HEAL (2026-08-27): re-ingest any
        # of the 2 days before that with ZERO ingested states — catches late
        # URMA availability and skipped Actions runs (the 2026-08-26
        # nationwide miss). Quiet national days are re-checked harmlessly.
        # Load discipline (2026-08-27): a day is (re)ingested ONLY while the DB
        # has zero states for it — the afternoon pass costs three count queries
        # unless something actually failed. No duplicate daily rewrites.
        # 2026-10-07: an unreadable count is an error (it used to skip the day).
        # 2026-10-07: counts only feed rows (obs_only=false): wind_obs_rescue's
        # 12:00Z rows must not make a deferred day look ingested. Days that
        # have rows are checked for SHORT WINDOWS (hour-aware completeness).
        today = dt.datetime.now(UTC).date()
        states_all = parse_states()
        for back in (1, 2, 3):
            d = today - dt.timedelta(days=back)
            try:
                n = run.rest_count("wind_days", {"valid_date": f"eq.{d.isoformat()}",
                                                 "obs_only": "is.false", "select": "state"})
            except Exception as e:
                run.error(d.isoformat(), "self-heal check", str(e), quarantine=False)
                continue
            if n == 0:
                if back > 1:
                    print(f"[self-heal] {d} has zero ingested states — re-running that date")
                dates.append(d.strftime("%Y%m%d"))
                continue
            try:
                short = stored_window_check(run, d.isoformat(), d.strftime("%Y%m%d"), states_all, step)
            except Exception as e:
                run.error(d.isoformat(), "completeness check", str(e), quarantine=False)
                continue
            if short and heal_from and d >= heal_from:
                print(f"[self-heal] {d} was written from a short URMA window — re-running it")
                dates.append(d.strftime("%Y%m%d"))
            elif short:
                run.warn(d.isoformat(), "STORED SHORT WINDOW: this day was written before its last "
                         "URMA hour(s) were published; not re-written (past-day repair needs owner "
                         "sign-off; HEAL_SHORT_WINDOWS_FROM enables it)")
        if not dates:
            print("All recent dates already ingested — nothing to do.")
    states = parse_states()

    print(f"URMA wind ingest ({conv}): {len(dates)} date(s), {len(states)} state(s), hour step {step}, "
          f"hours policy {policy}{' [BG_ONLY]' if bg_only else ''}{f' {flags}' if flags else ''}")
    total = 0
    for d in dates:
        if conv == "v4":
            total += process_date_v4(run, d, states, step, policy, flags, bg_only)
        else:
            total += process_date(run, d, states, step, policy, bg_only)

    if not run.dry_run and not bg_only:
        try:
            run._post("purge_old_wind", {"p_secret": run.secret}, 120)
            print("Purged wind data older than 2 years.")
        except Exception as e:
            print(f"[warn] purge failed: {e}")
    print(f"Done. {total} state-day(s) written across {len(dates)} date(s).")


if __name__ == "__main__":
    try:
        _run = fg.Run("wind")
    except fg.FeedError as e:
        print(f"::error::{e}")
        raise SystemExit(2)
    fg.main_guard(_run, main)
