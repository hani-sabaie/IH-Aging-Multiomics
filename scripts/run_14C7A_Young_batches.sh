#!/usr/bin/env bash
set -euo pipefail

RSCRIPT="/c/Program Files/R/R-4.5.2/bin/Rscript.exe"
SCRIPT="scripts/14C7A_CellChat_Full18_Batch_Audit.R"

OUTDIR="processed_results/12_CellChat/multiple_testing_audit/full18_nboot1000_BH/batches"
LOGDIR="$OUTDIR/logs"

mkdir -p "$LOGDIR"

TOTAL_LR=488
BATCH_SIZE=25
MAX_PARALLEL=2

run_batch() {
    local start="$1"
    local end="$2"

    local tag
    tag=$(printf "Young_LR%04d-%04d" "$start" "$end")

    local tests="$OUTDIR/${tag}_tests.tsv"
    local validation="$OUTDIR/${tag}_validation.tsv"
    local log="$LOGDIR/${tag}.log"

    if [[ -s "$tests" && -s "$validation" ]]; then
        echo "SKIP: $tag already complete"
        return 0
    fi

    echo "START: $tag"

    "$RSCRIPT" --vanilla \
        "$SCRIPT" \
        Young "$start" "$end" \
        >"$log" 2>&1

    if [[ ! -s "$tests" || ! -s "$validation" ]]; then
        echo "FAIL: expected output missing for $tag"
        return 1
    fi

    if ! grep -q $'\tTRUE$' "$validation"; then
        echo "FAIL: probability validation did not pass for $tag"
        tail -n 30 "$log"
        return 1
    fi

    echo "DONE: $tag"
}

running=0

for ((start=1; start<=TOTAL_LR; start+=BATCH_SIZE)); do

    end=$((start + BATCH_SIZE - 1))

    if (( end > TOTAL_LR )); then
        end=$TOTAL_LR
    fi

    run_batch "$start" "$end" &

    running=$((running + 1))

    if (( running >= MAX_PARALLEL )); then
        wait
        running=0
    fi
done

wait

echo
echo "============================================================"
echo "ALL YOUNG BATCHES COMPLETE"
echo "============================================================"
