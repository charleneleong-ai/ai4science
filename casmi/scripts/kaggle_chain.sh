#!/bin/bash
# Push a kernel snapshot as version V, wait for its commit run, gate it, submit it, and wait for the public score.
# Usage: KERNEL=<owner/slug> GATE=<gate.py> OUTDIR=<name> kaggle_chain.sh V SNAPSHOT_DIR AFTER_UTC_DATE MESSAGE
# Retries the push while the weekly GPU quota is spent (every 3 h) or both GPU sessions are busy (every 20 min).
# Only the log, CSV and JSON outputs are downloaded: the full output is large and the gate needs only these.
V=$1; NB=$2; AFTER=$3; MSG=$4
K=${KERNEL:?KERNEL is required}; C=enveda-CASMI26-molecule-id-mass-spectra
O=$HOME/casmi-data/kernel_outputs; D=$O/${OUTDIR:-v$V}; mkdir -p "$D"
k() { perl -e 'alarm 300; exec @ARGV' uvx kaggle "$@"; }

fails=0
while true; do
  p=$(k kernels push -p "$NB" 2>&1 | tail -1); echo "$p" | grep -q "successfully pushed" && break
  if echo "$p" | grep -q "GPU quota"; then echo "[$(date -u +%FT%TZ)] GPU quota exhausted; retrying in 3 h"; sleep 10800; continue; fi
  if echo "$p" | grep -q "session count"; then echo "[$(date -u +%FT%TZ)] GPU sessions busy; retrying in 20 min"; sleep 1200; continue; fi
  fails=$((fails+1)); [ $fails -ge 3 ] && break; sleep 60
done
echo "$p"
P=$(echo "$p" | grep -oE 'version [0-9]+' | grep -oE '[0-9]+')
[ "$P" = "$V" ] || { echo "pushed version '$P' is not $V - not submitting"; exit 1; }

sleep 90
until s=$(k kernels status "$K" 2>&1); echo "$s" | grep -qE 'COMPLETE|ERROR|CANCEL'; do sleep 120; done; echo "$s"
k kernels output "$K" -p "$D" --file-pattern '\.(log|csv|json)$' >/dev/null 2>&1
echo "$s" | grep -q COMPLETE || {
  grep -ohE '"data": *"[^"]*(Error|Traceback|SystemExit)[^"]*' "$D"/*.log | tail -5
  echo "RUN NOT COMPLETE - not submitting"; exit 1; }
uv run --no-project python "$HOME/casmi-data/gates/$GATE" "$D" || { echo "GATE FAIL"; exit 1; }

[ -n "$AFTER" ] && until [ "$(date -u +%Y-%m-%d)" \> "$AFTER" ]; do sleep 300; done
(cd "$D" && k competitions submit "$C" -f submission.csv -k "$K" -v "$V" -m "$MSG" 2>&1 | tail -1)
sleep 30; ref=""
for _ in 1 2 3 4 5; do
  ref=$(k competitions submissions "$C" --csv 2>/dev/null | grep -F "$MSG" | head -1 | cut -d, -f1)
  [ -n "$ref" ] && break; sleep 60
done
echo "ref $ref"
# hidden-test scoring has taken up to ~5 h; poll for 8 h
for _ in $(seq 1 96); do
  r=$(k competitions submissions "$C" --csv 2>/dev/null | grep "^$ref,")
  echo "$r" | grep -qE 'COMPLETE|ERROR' && break; sleep 300
done
echo "SCORE: $r"
