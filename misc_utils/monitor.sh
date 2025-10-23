#!/usr/bin/env bash

## You should not run this with /bin/time -v
## It will monitor the time command

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <command> [args...]" >&2
    exit 1
fi

OUTFILE="monitor.log"
> "$OUTFILE"

start_time=$(date +%s)

"$@" &
pid=$!
pgid=$(ps -o pgid= "$pid" | tr -d ' ')

cleanup() {
    echo
    kill -TERM -"$pgid" 2>/dev/null || true
    wait "$pid" 2>/dev/null || true
}
trap cleanup INT

while ps -p "$pid" > /dev/null 2>&1; do
    timestamp=$(date +%s.%N)
    mem_mb=$(ps -p "$pid" -o size= | awk '{printf "%.2f", $1/1024}')
    echo "$mem_mb" >> "$OUTFILE"
    sleep 0.1
done

# wait for process to finish to capture exit status
wait "$pid"
exit_code=$?

end_time=$(date +%s)
runtime=$(( end_time - start_time ))

echo "--------------------------------------------"
echo "Process finished at $(date +"%Y-%m-%d %H:%M:%S")"
echo "Exit code: $exit_code, Runtime: ${runtime}s"

python3 - <<'PYCODE'
import sys
import matplotlib.pyplot as plt

filename = 'monitor.log'
output = 'ram_stats.png'
interval = 0.1

times = []
mems = []

with open(filename, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if not parts:
            continue
        if len(parts) == 1:
            mems.append(float(parts[0]))
            times.append(i * interval)
        else:
            times.append(float(parts[0]))
            mems.append(float(parts[1]))

plt.figure(figsize=(8, 4))
plt.plot(times, mems, color='tab:blue', linewidth=1)
plt.title("RAM Usage Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Memory (MB)")
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(output, bbox_inches='tight')
PYCODE

rm -f "$OUTFILE"
