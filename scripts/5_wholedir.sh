#!/bin/bash
# Usage: 5_wholedir.sh SES [SITE]
#
# Run the pipeline over every subject in dicoms/uncompressed. The site and the
# subject label come from the DICOM headers via 0_identify_dicoms.py, not from
# the folder name -- folder names are inconsistent across sites and sometimes
# wrong (the folder "1-845" holds UCL participant P00<nnnn> -> rl3399).
#
# SITE is an optional filter: pass it to convert one site's subjects only.
#
# Skips subjects that already have a defaced T1w, and skips the second and
# subsequent folders of any duplicated study.

SES=$1
SITE_FILTER=$2

if [ -z "$SES" ]; then
    echo "usage: $(basename "$0") SES [SITE]   (SITE: UCB|UCL|SUT)" >&2
    exit 1
fi

ROOT=${BD2_ROOT:-/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged}
SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MANIFEST=$(mktemp -t bd2_manifest)
trap 'rm -f "$MANIFEST"' EXIT

echo "Identifying DICOM folders..."
python "$SCRIPTS/0_identify_dicoms.py" --tsv "$MANIFEST" >/dev/null || exit 1

echo "Session: $SES${SITE_FILTER:+  Site filter: $SITE_FILTER}"

SEEN_UIDS=""

# Fields are tab-separated and folder names may contain spaces, so read with
# IFS set to tab only. tail -n +2 drops the header row.
#
# The study UID field is read into STUDY_UID, not UID: bash makes UID readonly,
# so `read` into it fails at runtime with "readonly variable" and the loop
# silently processes nothing.
while IFS=$'\t' read -r FOLDER SITE SUB DATE N_SERIES N_FILES PID PNAME SCANNER SOFTWARE STUDY_UID FLAGS; do
    [ -z "$FOLDER" ] && continue

    if [ "$SITE" = "?" ] || [ "$SUB" = "?" ]; then
        echo "Skipping $FOLDER: could not identify (${FLAGS:-no detail})"
        continue
    fi

    if [ -n "$SITE_FILTER" ] && [ "$SITE" != "$SITE_FILTER" ]; then
        continue
    fi

    # One conversion per study, even if the same scan sits in two folders.
    case " $SEEN_UIDS " in
        *" $STUDY_UID "*)
            echo "Skipping $FOLDER: duplicate of an already-processed study ($SUB)"
            continue ;;
    esac
    SEEN_UIDS="$SEEN_UIDS $STUDY_UID"

    # Skip only if a defaced T1w is actually present. Testing for anat/ instead
    # -- as this did originally -- skips subjects that died mid-pipeline,
    # because 2_nifti_to_bids_naming.py creates the empty subfolders up front.
    shopt -s nullglob
    DONE=("$ROOT/bids/sub-$SUB/$SES/anat/"*_T1w.nii.gz)
    shopt -u nullglob
    if [ ${#DONE[@]} -gt 0 ]; then
        echo "Skipping sub-$SUB ($SITE): already has ${#DONE[@]} T1w"
        continue
    fi

    echo
    echo "=== Processing sub-$SUB  site $SITE  from $FOLDER ==="
    bash "$SCRIPTS/4_single_sub.sh" "$SITE" "$SUB" "$SES" \
         --src "$ROOT/dicoms/uncompressed/$FOLDER"
done < <(tail -n +2 "$MANIFEST")

# Refresh participants.tsv so the site column covers whatever just converted.
echo
python "$SCRIPTS/0_identify_dicoms.py" --participants >/dev/null
