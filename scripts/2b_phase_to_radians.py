import glob
import json
import os
import sys

import nibabel as nib
import numpy as np

# Rescale MEGRE phase images from scanner units to radians.
#
# Usage: python 2b_phase_to_radians.py SUB SES [--dry-run]
#
# dcm2niix writes Siemens phase in raw 12-bit scanner units (-4096..4095), not
# radians. Every QSM pipeline -- QSMxT, SEPIA, TGV-QSM, MEDI -- expects
# -pi..+pi. Handed scanner units they do not error; they unwrap "successfully"
# and produce a susceptibility map that is wrong by ~1300x. Nothing downstream
# catches this except sanity-checking the final chi values, so it gets fixed
# here, once, right after renaming.
#
# Runs after 2_nifti_to_bids_naming.py and before 3_deface.sh, but it searches
# recursively so it works either side of the sort into anat/.

BIDS_DIR = os.environ.get(
    "BD2_ROOT", "/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged"
) + "/bids"

# Anything with |phase| below this is already in radians. pi plus a little
# slack for interpolation overshoot.
RADIAN_CEILING = 3.25

# Nominal ranges, keyed by whether the stored data is signed. These are the
# protocol values, not the observed extremes -- a volume whose phase never
# quite reaches the rail would otherwise be rescaled by a slightly wrong
# factor, and that error is invisible in the output.
SIGNED_HALF_RANGE = 4096.0     # -4096..4095  -> center 0
UNSIGNED_HALF_RANGE = 2048.0   #     0..4095  -> center 2048


def classify(lo, hi):
    """Return (center, half_range, label) for observed phase extrema.

    Raises ValueError on anything unrecognised. Refusing loudly beats guessing:
    a wrong scale factor here is silent all the way to the chi map.
    """
    if -RADIAN_CEILING <= lo and hi <= RADIAN_CEILING:
        return None, None, "radians"
    if lo < -100.0:
        if lo < -1.05 * SIGNED_HALF_RANGE or hi > 1.05 * SIGNED_HALF_RANGE:
            raise ValueError("range [%g, %g] exceeds signed 12-bit" % (lo, hi))
        return 0.0, SIGNED_HALF_RANGE, "signed 12-bit"
    if hi > 100.0:
        if hi > 1.05 * 2 * UNSIGNED_HALF_RANGE:
            raise ValueError("range [%g, %g] exceeds unsigned 12-bit" % (lo, hi))
        return UNSIGNED_HALF_RANGE, UNSIGNED_HALF_RANGE, "unsigned 12-bit"
    raise ValueError("range [%g, %g] matches no known phase encoding" % (lo, hi))


def sidecar_says_radians(nii_path):
    """True if this image's JSON already records radians.

    Cheaper than loading the volume, and it makes the script idempotent even
    for a subject whose phase happens to be small.
    """
    json_path = nii_path[: -len(".nii.gz")] + ".json"
    if not os.path.exists(json_path):
        return False
    with open(json_path) as fh:
        return json.load(fh).get("Units") == "rad"


def mark_sidecar(nii_path, dry_run):
    """Record the units in the JSON so the conversion is self-documenting."""
    json_path = nii_path[: -len(".nii.gz")] + ".json"
    if not os.path.exists(json_path):
        print("    WARNING: no sidecar at %s" % os.path.basename(json_path))
        return
    if dry_run:
        return
    with open(json_path) as fh:
        meta = json.load(fh)
    meta["Units"] = "rad"
    with open(json_path, "w") as fh:
        json.dump(meta, fh, indent=4)
        fh.write("\n")


def rescale(nii_path, dry_run):
    """Convert one phase image to radians in place. Returns a status string."""
    name = os.path.basename(nii_path)

    if sidecar_says_radians(nii_path):
        print("  skip    %s (sidecar already says rad)" % name)
        return "skipped"

    img = nib.load(nii_path)
    data = img.get_fdata()
    lo, hi = float(data.min()), float(data.max())

    try:
        center, half_range, label = classify(lo, hi)
    except ValueError as exc:
        print("  REFUSE  %s: %s" % (name, exc))
        return "refused"

    if label == "radians":
        print("  skip    %s (already radians, [%.3f, %.3f])" % (name, lo, hi))
        mark_sidecar(nii_path, dry_run)
        return "skipped"

    print("  rescale %s  %s [%.0f, %.0f] -> rad" % (name, label, lo, hi))
    if dry_run:
        return "rescaled"

    radians = (data - center) * (np.pi / half_range)

    # Build the output on a fresh header. Inheriting the input's header would
    # carry over its scl_slope/scl_inter (2 and -4096 for Siemens int16), and
    # those would be reapplied on the next read.
    out = nib.Nifti1Image(radians.astype(np.float32), img.affine)
    out.header.set_xyzt_units(*img.header.get_xyzt_units())
    out.set_data_dtype(np.float32)
    out.to_filename(nii_path)

    mark_sidecar(nii_path, dry_run)
    return "rescaled"


def main():
    args = [a for a in sys.argv[1:] if a != "--dry-run"]
    dry_run = "--dry-run" in sys.argv

    if len(args) != 2:
        sys.exit("usage: %s SUB SES [--dry-run]" % os.path.basename(sys.argv[0]))

    sub, ses = args
    session = os.path.join(BIDS_DIR, "sub-" + sub, ses)
    if not os.path.isdir(session):
        sys.exit("no such session: %s" % session)

    pattern = os.path.join(session, "**", "*part-phase*MEGRE.nii.gz")
    phases = sorted(glob.glob(pattern, recursive=True))

    print("phase -> radians for sub-%s %s" % (sub, ses))
    if not phases:
        print("  no MEGRE phase images found -- nothing to do")
        return

    counts = {"rescaled": 0, "skipped": 0, "refused": 0}
    for path in phases:
        counts[rescale(path, dry_run)] += 1

    print("  %d rescaled, %d skipped, %d refused%s"
          % (counts["rescaled"], counts["skipped"], counts["refused"],
             "  (dry run -- nothing written)" if dry_run else ""))

    if counts["refused"]:
        sys.exit(1)


main()
