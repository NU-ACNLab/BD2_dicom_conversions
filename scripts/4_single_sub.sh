#!/bin/bash
# Usage: 4_single_sub.sh SITE SUB SES [--dry-run] [--src DIR]
#
# Full pipeline for one subject/session. --dry-run is passed to the renamer and
# the phase rescaler; conversion still runs and QC is skipped.
#
# --src names the source DICOM folder when it is not called SUB -- see
# 1_dic_to_nifti.sh. 5_wholedir.sh fills it in from the DICOM headers.

SITE=$1
SUB=$2
SES=$3
shift 3 2>/dev/null

DRY_RUN=""
SRC=""
while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN="--dry-run"; shift ;;
        --src) SRC=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

if [ -z "$SITE" ] || [ -z "$SUB" ] || [ -z "$SES" ]; then
    echo "usage: $(basename "$0") SITE SUB SES [--dry-run] [--src DIR]   (SITE: UCB|UCL|SUT)" >&2
    exit 1
fi

SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SRC_ARGS=()
[ -n "$SRC" ] && SRC_ARGS=(--src "$SRC")

bash "$SCRIPTS/1_dic_to_nifti.sh" "$SITE" "$SUB" "$SES" "${SRC_ARGS[@]}" || exit 1
python "$SCRIPTS/2_nifti_to_bids_naming.py" "$SITE" "$SUB" "$SES" $DRY_RUN || exit 1
# Must run after the renamer -- it identifies phase by the part-phase entity,
# which only exists once 2_ has run. Safe either side of 3_'s sort into anat/.
python "$SCRIPTS/2b_phase_to_radians.py" "$SUB" "$SES" $DRY_RUN || exit 1
bash "$SCRIPTS/3_deface.sh" "$SITE" "$SUB" "$SES" || exit 1

# QC is advisory here: report problems but do not fail the subject, so a batch
# run gets through. Both scripts exit non-zero on trouble if you want to gate.
# Skipped under --dry-run, where there is nothing renamed to inspect.
if [ -z "$DRY_RUN" ]; then
    python "$SCRIPTS/6_qc_qsm.py" "$SES" "$SUB" || true
    python "$SCRIPTS/7_qc_images.py" "$SES" "$SUB" || true
fi
