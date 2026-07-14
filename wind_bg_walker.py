#!/usr/bin/env python3
"""
Backfills EXACT URMA station backgrounds + coarse field samples for the last
2 years (BG_ONLY mode: downloads each date's 24 hourly GUST slices, samples
them, writes hz_station_bg + hz_bg_coarse, skips the swath rebuild). Walks
backward from yesterday, cursor-saved per date (key 'anlbg'), so it resumes
across runs and never repeats work. ~2 min/date; full depth in ~2 days of
scheduled runs.
Env: SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: TIME_BUDGET_MIN (100), START_DATE/END_DATE (YYYYMMDD)
"""
import os, json, time, datetime as dt
import requests
os.environ["BG_ONLY"] = "1"
import wind_ingest as w


def rpc(base, anon, name, payload):
    r = requests.post(f"{base}/rest/v1/rpc/{name}",
        headers={"apikey": anon, "Authorization": f"Bearer {anon}",
                 "Content-Type": "application/json"},
        data=json.dumps(payload), timeout=60)
    if r.status_code >= 300:
        raise RuntimeError(f"{name} {r.status_code}: {r.text[:200]}")
    return r.json() if r.text and r.text != "null" else None


def main():
    t0 = time.time()
    budget = 60 * int(os.environ.get("TIME_BUDGET_MIN", "100"))
    base = os.environ["SUPABASE_URL"].rstrip("/")
    anon = os.environ["SUPABASE_ANON_KEY"]
    secret = os.environ["INGEST_SECRET"]
    today = dt.date.today()
    end = dt.datetime.strptime(os.environ["END_DATE"], "%Y%m%d").date() \
        if os.environ.get("END_DATE") else today - dt.timedelta(days=1)
    start = dt.datetime.strptime(os.environ["START_DATE"], "%Y%m%d").date() \
        if os.environ.get("START_DATE") else today - dt.timedelta(days=730)
    cur = rpc(base, anon, "hz_backfill_get", {"p_key": "anlbg"})
    cursor = dt.datetime.strptime(cur, "%Y-%m-%d").date() if cur \
        else end + dt.timedelta(days=1)
    print(f"ANL-bg walker. Budget {budget//60} min. Resuming before {cursor}.")
    while True:
        day = cursor - dt.timedelta(days=1)
        if day < start:
            print(f"ANL background backfill COMPLETE: reached {start}."); break
        if time.time() - t0 > budget:
            print(f"Budget reached. Next run resumes before {cursor}."); break
        try:
            w.process_date(day.strftime("%Y%m%d"), [], base, anon, secret, 1)
            print(f"  {day}: backgrounds written [{int(time.time()-t0)}s]")
        except Exception as e:
            print(f"  [error] {day}: {e} -- advancing past it")
        rpc(base, anon, "hz_backfill_set",
            {"p_secret": secret, "p_key": "anlbg",
             "p_value": day.strftime("%Y-%m-%d")})
        cursor = day


if __name__ == "__main__":
    main()
