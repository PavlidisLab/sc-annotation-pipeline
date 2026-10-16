#!/usr/bin/env bash
# Runs the pipeline's test profiles sequentially and reports pass/fail for each.
# Usage:
#   scripts/run_test_profiles.sh                  # run all 8 profiles
#   scripts/run_test_profiles.sh test_mouse test_human   # run only these

set -uo pipefail

ALL_PROFILES=(
    test_mouse
    test_human
    test_mouse_persample
    test_human_persample
    test_mouse_local
    test_human_local
    test_mouse_upload_off
    test_human_upload_off
)

PROFILES=("$@")
if [ ${#PROFILES[@]} -eq 0 ]; then
    PROFILES=("${ALL_PROFILES[@]}")
fi

declare -A RESULTS

for profile in "${PROFILES[@]}"; do
    echo "=== Running profile: ${profile} ==="
    if nextflow run main.nf -profile "${profile}",singularity --upload_gemma false -work-dir "work/test_${profile}" -resume; then
        RESULTS["${profile}"]="PASS"
    else
        RESULTS["${profile}"]="FAIL"
    fi
    echo
done

echo "=== Summary ==="
overall_status=0
for profile in "${PROFILES[@]}"; do
    printf '%-25s %s\n' "${profile}" "${RESULTS[${profile}]}"
    if [ "${RESULTS[${profile}]}" != "PASS" ]; then
        overall_status=1
    fi
done

if [ "${overall_status}" -eq 0 ]; then
    echo "All test profiles passed - safe to tag a release."
else
    echo "One or more test profiles failed - do not tag a release."
fi

exit "${overall_status}"
