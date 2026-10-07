#!/usr/bin/env python3
"""
Timeout-proof self-walking backfill for the wind tool.

Design: instead of fixed-size chunks (which can run long on hurricane weeks),
each run processes dates ONE AT A TIME against a wall-clock budget
(TIME_BUDGET_MIN, default 100). Progress is saved to Supabase after EVERY
COMPLETED date, so a run can never hit the workflow timeout with unsaved work,
and the next scheduled run resumes exactly where the last stopped. It walks
BACKWARD from yesterday to 2 years ago (newest history fills first) at effort 1
(full hourly resolution), all permitted states, HOURS_POLICY=strict.

2026-10-07 (Archive Phase 5, Stage 1): the walker STOPS on the first date that
did not ingest completely and does NOT move the cursor past it (it used to log
and advance, leaving permanent silent holes). The run exits non-zero.
This workflow is disabled (Connor); this only matters if it is re-enabled.

Env: SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: TIME_BUDGET_MIN (default 100), START_DATE/END_DATE (YYYYMMDD), DRY_RUN
"""
import os, time, datetime as dt
import feedguard as fg
import wind_ingest as w   # reuse the tested ingester


def walk(run):
    t0 = time.time()
    budget_s = 60 * int((os.environ.get("TIME_BUDGET_MIN") or "100").strip() or "100")

    today = dt.date.today()
    end   = dt.datetime.strptime(os.environ["END_DATE"], "%Y%m%d").date() \
            if os.environ.get("END_DATE") else today - dt.timedelta(days=1)
    start = dt.datetime.strptime(os.environ["START_DATE"], "%Y%m%d").date() \
            if os.environ.get("START_DATE") else today - dt.timedelta(days=730)

    cur = run.rpc_read("backfill_get", {"p_key": "wind", "p_secret": run.secret}).json()
    cursor = dt.datetime.strptime(cur, "%Y-%m-%d").date() if cur else end + dt.timedelta(days=1)

    states = sorted(w.PERMITTED_STATES)
    done = 0
    print(f"Walker start. Budget {budget_s//60} min. Resuming before {cursor}. Target floor {start}.")

    while True:
        day = cursor - dt.timedelta(days=1)
        if day < start:
            print(f"Backfill COMPLETE: reached {start}.")
            break
        elapsed = time.time() - t0
        if elapsed > budget_s:
            print(f"Time budget reached after {done} date(s). Next run resumes before {cursor}.")
            break
        n = w.process_date(run, day.strftime("%Y%m%d"), states, 1, "strict")
        if run.day(day.isoformat())["status"] == "error":
            print(f"  [STOP] {day} did not ingest completely; cursor stays at {cursor}.")
            break
        print(f"  {day}: {n} state-day(s) written  [{int(elapsed)}s elapsed]")
        if not run.dry_run:
            run._post("backfill_set", {"p_key": "wind", "p_value": day.strftime("%Y-%m-%d"),
                                       "p_secret": run.secret}, 60)
        cursor = day
        done += 1


if __name__ == "__main__":
    try:
        _run = fg.Run("wind-backfill")
    except fg.FeedError as e:
        print(f"::error::{e}")
        raise SystemExit(2)
    fg.main_guard(_run, walk)
