#!/usr/bin/env python3
"""
Timeout-proof self-walking backfill for the wind tool.

Design: instead of fixed-size chunks (which can run long on hurricane weeks),
each run processes dates ONE AT A TIME against a wall-clock budget
(TIME_BUDGET_MIN, default 100). Progress is saved to Supabase after EVERY
date, so a run can never hit the workflow timeout with unsaved work, and the
next scheduled run resumes exactly where the last stopped. It walks BACKWARD
from yesterday to 2 years ago (newest history fills first) at effort 1
(full hourly resolution), all permitted states.

Env: SUPABASE_URL, SUPABASE_ANON_KEY, INGEST_SECRET
Optional: TIME_BUDGET_MIN (default 100), START_DATE/END_DATE (YYYYMMDD)
"""
import os, json, time, datetime as dt
import requests
import wind_ingest as w   # reuse the tested ingester


def rpc(base, anon, name, payload):
    r = requests.post(f"{base}/rest/v1/rpc/{name}",
        headers={"apikey": anon, "Authorization": f"Bearer {anon}",
                 "Content-Type": "application/json"},
        data=json.dumps(payload), timeout=60)
    if r.status_code >= 300:
        raise RuntimeError(f"{name} {r.status_code}: {r.text[:200]}")
    return r.json() if r.text else None


def main():
    t0 = time.time()
    budget_s = 60 * int((os.environ.get("TIME_BUDGET_MIN") or "100").strip() or "100")

    base   = os.environ["SUPABASE_URL"].rstrip("/")
    anon   = os.environ["SUPABASE_ANON_KEY"]
    secret = os.environ["INGEST_SECRET"]

    today = dt.date.today()
    end   = dt.datetime.strptime(os.environ["END_DATE"], "%Y%m%d").date() \
            if os.environ.get("END_DATE") else today - dt.timedelta(days=1)
    start = dt.datetime.strptime(os.environ["START_DATE"], "%Y%m%d").date() \
            if os.environ.get("START_DATE") else today - dt.timedelta(days=730)

    cur = rpc(base, anon, "backfill_get", {"p_key": "wind"})
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
        ds = day.strftime("%Y%m%d")
        try:
            n = w.process_date(ds, states, base, anon, secret, step=1)
            print(f"  {day}: {n} state-day(s) written  [{int(elapsed)}s elapsed]")
        except Exception as e:
            # log and still advance: a permanently-missing archive day must not
            # wedge the walker. Transient errors are already retried inside
            # process_date's download/rpc layers.
            print(f"  [error] {day}: {e} -- advancing past it")
        rpc(base, anon, "backfill_set", {"p_key": "wind", "p_value": day.strftime("%Y-%m-%d")})
        cursor = day
        done += 1


if __name__ == "__main__":
    main()
