# stormauditor-wind-feed

## Operations (Archive Phase 5 Stage 1, 2026-10-07)
- `DATE=YYYY-MM-DD` or `DATE=START..END` (also `YYYYMMDD`, `a:b`, comma lists; workflow input `date`). Blank = scheduled mode.
- **Hours rule.** A day is written only when every UTC hour it needs is published on `noaa-urma-pds`: each state's local-day window (23/24/25 h on DST days) plus the 06Z-06Z background window. URMA lands about 7 h after valid time. Scheduled runs (`HOURS_POLICY=defer`) therefore write nothing for yesterday at 10:10Z and 12:10Z, and the 16:10Z dispatch writes it complete. Explicit dates use `strict` (a missing hour fails). `legacy` reproduces the old short windows and exists only for parity tests.
- **Self-heal.** It ignores `obs_only` rescue rows and reports stored days written before their last hour was published (short window). It never re-writes them unless `HEAL_SHORT_WINDOWS_FROM=YYYY-MM-DD` is set (owner sign-off).
- URMA before 2024-04-18 has no `.idx`. The gust message is found by GRIB metadata through ranged header reads, and the bytes are identical to the `.idx` path.
- `DRY_RUN=1`, `FEED_RESULT {json}`, quarantine and exit codes work as in the hail feed (`feedguard.py`, shared verbatim with the other feed repos).
