#!/bin/bash
# Usage: 3_deface.sh SITE SUB SES
#
# Sort the renamed files into BIDS subfolders, deface every T1w, and clear out
# the scans we don't keep. Only the cleanup patterns differ by site.

SITE=$1
SUB=$2
SES=$3

if [ -z "$SITE" ] || [ -z "$SUB" ] || [ -z "$SES" ]; then
    echo "usage: $(basename "$0") SITE SUB SES   (SITE: UCB|UCL|SUT)" >&2
    exit 1
fi

ROOT=${BD2_ROOT:-/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged}
BASE=$ROOT/bids/sub-$SUB/$SES

# Scans discarded after conversion. UCB removed localizers; SUT removed
# localizers and MPR reconstructions; UCL swept anything still prefixed with
# the raw subject ID.
case "$SITE" in
    UCB) CLEANUP=("*localizer*") ;;
    SUT) CLEANUP=("*localizer*" "*MPR*") ;;
    UCL) CLEANUP=("*localizer*" "${SUB}-*") ;;
    *) echo "unknown site: $SITE" >&2; exit 1 ;;
esac

echo "Site: $SITE  Subject: $SUB  Session: $SES"

# Sort into BIDS subfolders. Order matters: the QSM echoes carry
# part-mag/part-phase entities, so the *MEGRE* rule must claim them for anat/
# before the *phase* and *mag* rules would sweep them into fmap/.
mv "$BASE"/*T1w.* "$BASE"/anat/ 2>/dev/null
mv "$BASE"/*MEGRE* "$BASE"/anat/ 2>/dev/null
mv "$BASE"/*phase* "$BASE"/fmap/ 2>/dev/null
mv "$BASE"/*mag* "$BASE"/fmap/ 2>/dev/null
mv "$BASE"/*epi* "$BASE"/fmap/ 2>/dev/null
mv "$BASE"/*dwi* "$BASE"/dwi/ 2>/dev/null
mv "$BASE"/*bold* "$BASE"/func/ 2>/dev/null

# Deface every T1w present. The originals hardcoded run-1/run-2; SUT only ever
# defaced run-1, so a second T1w shipped with the face intact.
shopt -s nullglob
for T1W in "$BASE"/anat/*_T1w.nii.gz; do
    echo "Defacing $T1W"
    pydeface --verbose "$T1W"

    DEFACED=${T1W%.nii.gz}_defaced.nii.gz
    if [ -f "$DEFACED" ]; then
        mv "$DEFACED" "$T1W"
    else
        echo "WARNING: pydeface produced no output for $T1W -- left as is" >&2
    fi
done

for PATTERN in "${CLEANUP[@]}"; do
    rm -f "$BASE"/$PATTERN
done
shopt -u nullglob
