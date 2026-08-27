import argparse
import glob
import json
import os
import re
import sys

# Fieldmap naming + the sidecar metadata that BIDS requires and dcm2niix does
# not write. Safe to re-run after every conversion; it is idempotent.
#
# Usage: python 9_bids_metadata.py [SUB ...] [--dry-run]
#
#   python 9_bids_metadata.py             # whole dataset
#   python 9_bids_metadata.py rb860       # one subject
#   python 9_bids_metadata.py --dry-run   # report only
#
# Three things:
#
# 1. Renames legacy fieldmaps. `task` is not a valid entity for fmap/ -- the
#    BIDS filename template for a pepolar epi is
#      sub-<>[_ses-<>][_acq-<>][_ce-<>]_dir-<>[_run-<>][_part-<>][_chunk-<>]_epi
#    so `task-mid_run-01_dir-AP_epi` was invalid twice over: a `task` entity
#    that does not belong, and `run` before `dir` when the order is fixed.
#    Becomes `acq-mid_dir-AP_epi`. `run` is dropped rather than reordered
#    because there is exactly one fieldmap pair per task. dwi carried the same
#    entity-order bug (`run-1_dir-AP_dwi`) and gets the same treatment.
#
# 2. Writes IntendedFor into every fieldmap sidecar, as paths relative to the
#    subject directory (the long-standing form, resolved by every fMRIPrep
#    version). Without it nothing applies susceptibility distortion correction.
#
# 3. Writes TaskName into every func sidecar. BIDS requires it and the
#    validator errors without it. The value is the task entity from the
#    filename verbatim -- the validator compares the two, so a prettier
#    "Monetary Incentive Delay" would trip a mismatch warning.

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIDS_DIR = os.path.join(ROOT, "bids")

# task entity in the bold filenames <-> acq entity in the fieldmap filenames
TASKS = ("mid", "RL")

LEGACY_FMAP = re.compile(
    r"^(?P<sub>sub-[^_]+)_(?P<ses>ses-[^_]+)_task-(?P<task>[^_]+)_run-\d+_"
    r"dir-(?P<dir>[^_]+)_epi(?P<ext>\.nii\.gz|\.json)$")

# dwi had the same entity-order bug: run before dir, when the template fixes
# dir first. Same fix, and run is dropped for the same reason.
LEGACY_DWI = re.compile(
    r"^(?P<sub>sub-[^_]+)_(?P<ses>ses-[^_]+)_run-\d+_"
    r"dir-(?P<dir>[^_]+)_dwi(?P<ext>\.nii\.gz|\.json|\.bval|\.bvec)$")

EXTS = (".nii.gz", ".json")


def subjects(requested):
    found = sorted(os.path.basename(p) for p in glob.glob(BIDS_DIR + "/sub-*")
                   if os.path.isdir(p))
    if not requested:
        return found
    labels = {s: "sub-" + s.replace("sub-", "") for s in requested}
    missing = [s for s, full in labels.items() if full not in found]
    if missing:
        sys.exit("no such subject in bids/: " + ", ".join(missing))
    return [labels[s] for s in requested]


def sessions(sub):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(BIDS_DIR, sub, "ses-*"))
                  if os.path.isdir(p))


def rename_legacy(sub, ses, modality, pattern, template, dry_run, log):
    """Rename files matching a legacy pattern to the BIDS-valid form."""
    folder = os.path.join(BIDS_DIR, sub, ses, modality)
    if not os.path.isdir(folder):
        return
    for name in sorted(os.listdir(folder)):
        m = pattern.match(name)
        if not m:
            continue
        new = template % m.groupdict()
        src, dst = os.path.join(folder, name), os.path.join(folder, new)
        if os.path.exists(dst):
            log.append(("SKIP", "%s/%s/%s/%s -> %s already exists"
                        % (sub, ses, modality, name, new)))
            continue
        if not dry_run:
            os.rename(src, dst)
        log.append(("rename", "%s/%s/%s/%s -> %s" % (sub, ses, modality, name, new)))


def rename_fmaps(sub, ses, dry_run, log):
    """task-<t>_run-01_dir-<d>_epi -> acq-<t>_dir-<d>_epi, nifti and sidecar."""
    rename_legacy(sub, ses, "fmap", LEGACY_FMAP,
                  "%(sub)s_%(ses)s_acq-%(task)s_dir-%(dir)s_epi%(ext)s", dry_run, log)


def rename_dwi(sub, ses, dry_run, log):
    """run-1_dir-<d>_dwi -> dir-<d>_dwi, nifti, sidecar, bval and bvec."""
    rename_legacy(sub, ses, "dwi", LEGACY_DWI,
                  "%(sub)s_%(ses)s_dir-%(dir)s_dwi%(ext)s", dry_run, log)


def bolds_for_task(sub, ses, task):
    """Subject-relative paths of every bold run of one task, sorted."""
    pattern = os.path.join(BIDS_DIR, sub, ses, "func",
                           "%s_%s_task-%s_*_bold.nii.gz" % (sub, ses, task))
    return sorted("%s/func/%s" % (ses, os.path.basename(p)) for p in glob.glob(pattern))


def update_json(path, updates, dry_run, log, label):
    """Merge keys into a sidecar, only writing when something actually changes."""
    with open(path) as fo:
        data = json.load(fo)
    changed = {k: v for k, v in updates.items() if data.get(k) != v}
    if not changed:
        return False
    data.update(changed)
    if not dry_run:
        with open(path, "w") as fo:
            json.dump(data, fo, indent=4)
            fo.write("\n")
    log.append((label, "%s  %s" % (os.path.relpath(path, BIDS_DIR),
                                   ", ".join(sorted(changed)))))
    return True


def write_intended_for(sub, ses, dry_run, log):
    fmap = os.path.join(BIDS_DIR, sub, ses, "fmap")
    if not os.path.isdir(fmap):
        return
    for task in TASKS:
        targets = bolds_for_task(sub, ses, task)
        sidecars = sorted(glob.glob(os.path.join(
            fmap, "%s_%s_acq-%s_dir-*_epi.json" % (sub, ses, task))))
        if not sidecars:
            continue
        if not targets:
            log.append(("WARN", "%s/%s fmap acq-%s has no task-%s bold to point at"
                        % (sub, ses, task, task)))
            continue
        for path in sidecars:
            update_json(path, {"IntendedFor": targets}, dry_run, log, "IntendedFor")


def write_task_name(sub, ses, dry_run, log):
    func = os.path.join(BIDS_DIR, sub, ses, "func")
    if not os.path.isdir(func):
        return
    for path in sorted(glob.glob(os.path.join(func, "*_bold.json"))):
        m = re.search(r"_task-([^_]+)_", os.path.basename(path))
        if not m:
            log.append(("WARN", "no task entity in " + os.path.basename(path)))
            continue
        update_json(path, {"TaskName": m.group(1)}, dry_run, log, "TaskName")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("subjects", nargs="*", help="subject labels (default: all)")
    ap.add_argument("--dry-run", action="store_true", help="report only, change nothing")
    args = ap.parse_args()

    log = []
    for sub in subjects(args.subjects):
        for ses in sessions(sub):
            rename_fmaps(sub, ses, args.dry_run, log)
            rename_dwi(sub, ses, args.dry_run, log)
            write_intended_for(sub, ses, args.dry_run, log)
            write_task_name(sub, ses, args.dry_run, log)

    counts = {}
    for kind, _ in log:
        counts[kind] = counts.get(kind, 0) + 1
    for kind in ("rename", "IntendedFor", "TaskName", "SKIP", "WARN"):
        if kind in counts:
            print("%-12s %d" % (kind, counts[kind]))
    if not log:
        print("nothing to do -- dataset already up to date")
    for kind, msg in log:
        if kind in ("SKIP", "WARN"):
            print("  %s %s" % (kind, msg))
    if args.dry_run:
        print("DRY RUN -- nothing written")


if __name__ == "__main__":
    main()
