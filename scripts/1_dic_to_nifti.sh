#!/bin/bash
# Usage: 1_dic_to_nifti.sh SITE SUB SES [--src DIR]
#
# DICOM -> NIfTI. Identical for all three sites; SITE is accepted so every
# script in the pipeline takes the same arguments in the same order, and so
# per-site dcm2niix options have somewhere to go if that ever changes.
#
# --src names the source DICOM folder when it is not called SUB. Sites export
# under whatever name they like ("1-845", "rs1403_MR_1"), while SUB is the
# study label derived from the headers by 0_identify_dicoms.py. Defaults to
# dicoms/uncompressed/SUB, so a subject whose folder already matches its label
# needs no flag.

SITE=$1
SUB=$2
SES=$3
shift 3 2>/dev/null

SRC=""
while [ $# -gt 0 ]; do
    case "$1" in
        --src) SRC=$2; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
done

if [ -z "$SITE" ] || [ -z "$SUB" ] || [ -z "$SES" ]; then
    echo "usage: $(basename "$0") SITE SUB SES [--src DIR]   (SITE: UCB|UCL|SUT)" >&2
    exit 1
fi

ROOT=${BD2_ROOT:-/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged}

SRC=${SRC:-$ROOT/dicoms/uncompressed/$SUB}

if [ ! -d "$SRC" ]; then
    echo "no such DICOM folder: $SRC" >&2
    echo "run 0_identify_dicoms.py to see folder -> subject mapping" >&2
    exit 1
fi

echo "Site: $SITE  Subject: $SUB  Session: $SES"
echo "Source: $SRC"

case "$SITE" in
    UCB|UCL|SUT) DCM2NIIX_OPTS=(-b y -z o -w 1 -f '%n--%d--s%s--e%e') ;;
    *) echo "unknown site: $SITE" >&2; exit 1 ;;
esac

OUTPUT=$ROOT/bids/sub-$SUB/$SES
mkdir -p "$OUTPUT"

# Convert each series independently. A series that dcm2niix cannot use must not
# abort the subject: Siemens exports a PhoenixZIPReport series (modality SR, a
# protocol dump with no images) that dcm2niix exits 2 on, and because it is
# numbered 99 it sorts last, so its failure used to become this script's exit
# code. 4_single_sub.sh then stopped before renaming and the subject was left
# as raw dcm2niix output.
FAILED=()
for SCAN in "$SRC"/*; do
    [ -d "$SCAN" ] || continue
    echo "$SCAN"
    if ! dcm2niix "${DCM2NIIX_OPTS[@]}" -o "$OUTPUT" "$SCAN"; then
        echo "  no images in $(basename "$SCAN") -- skipping" >&2
        FAILED+=("$(basename "$SCAN")")
    fi
done

shopt -s nullglob
PRODUCED=("$OUTPUT"/*.nii.gz)
shopt -u nullglob

echo "Converted ${#PRODUCED[@]} NIfTI file(s) into $OUTPUT"
[ ${#FAILED[@]} -gt 0 ] && echo "Skipped non-image series: ${FAILED[*]}"

# Only a subject that produced nothing at all is a real failure.
if [ ${#PRODUCED[@]} -eq 0 ]; then
    echo "no NIfTI produced for $SUB -- aborting" >&2
    exit 1
fi
