#!/usr/bin/env python3
"""
Backfills EXACT URMA station backgrounds + coarse field samples for the last
2 years (BG_ONLY mode: downloads each date's 24 hourly GUST slices, samples
them, writes hz_station_bg + hz_bg_coarse, skips the swath rebuild). Walks
backward from yesterday, cursor-saved per COMPLETED date (key 'anlbg'), so it
resumes across runs and never repeats work. ~2 min/date.

2026-10-07 (Archive Phase 5, Stage 1): stops on the first incomplete date and
never advances the cursor past it (it used to log and advance). HOURS_POLICY
is strict. DISABLED by Connor 2026-09-16; re-enable only on his direction.
Env: SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: TIME_BUDGET_MIN (100), START_DATE/END_DATE (YYYYMMDD), DRY_RUN
"""
import os, time, datetime as dt
import feedguard as fg
import wind_ingest as w


def walk(run):
    t0 = time.time()
    budget = 60 * int(os.environ.get("TIME_BUDGET_MIN", "100"))
    today = dt.date.today()
    end = dt.datetime.strptime(os.environ["END_DATE"], "%Y%m%d").date() \
        if os.environ.get("END_DATE") else today - dt.timedelta(days=1)
    start = dt.datetime.strptime(os.environ["START_DATE"], "%Y%m%d").date() \
        if os.environ.get("START_DATE") else today - dt.timedelta(days=730)
    r = run.rpc_read("hz_backfill_get", {"p_key": "anlbg", "p_secret": run.secret})
    cur = r.json() if r.text and r.text != "null" else None
    cursor = dt.datetime.strptime(cur, "%Y-%m-%d").date() if cur \
        else end + dt.timedelta(days=1)
    print(f"ANL-bg walker. Budget {budget//60} min. Resuming before {cursor}.")
    while True:
        day = cursor - dt.timedelta(days=1)
        if day < start:
            print(f"ANL background backfill COMPLETE: reached {start}."); break
        if time.time() - t0 > budget:
            print(f"Budget reached. Next run resumes before {cursor}."); break
        w.process_date(run, day.strftime("%Y%m%d"), [], 1, "strict", bg_only=True)
        if run.day(day.isoformat())["status"] == "error":
            print(f"  [STOP] {day} did not ingest completely; cursor stays at {cursor}."); break
        print(f"  {day}: backgrounds written [{int(time.time()-t0)}s]")
        if not run.dry_run:
            run._post("hz_backfill_set", {"p_secret": run.secret, "p_key": "anlbg",
                                          "p_value": day.strftime("%Y-%m-%d")}, 60)
        cursor = day


if __name__ == "__main__":
    try:
        _run = fg.Run("anl-bg-backfill")
    except fg.FeedError as e:
        print(f"::error::{e}")
        raise SystemExit(2)
    fg.main_guard(_run, walk)
