import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from nilearn.glm.first_level import FirstLevelModel

# First-level GLM for the MID task, on fMRIPrep output.
#
# Usage: python 10_first_level_mid.py [SUB ...] [options]
#
#   python 10_first_level_mid.py                       # every subject
#   python 10_first_level_mid.py rs1356 rb1217         # just these
#   python 10_first_level_mid.py --dry-run             # report what it would fit
#   python 10_first_level_mid.py --fmriprep /path/...  # override paths
#
# Both MID runs are fit in ONE model rather than separately. There are only 6
# anticipation trials per condition per run, and as few as 1 Hit per run for
# rew_reward.neut, so per-run estimates would be very noisy; fitting jointly
# doubles the trials behind every contrast and yields one map per subject
# ready for a group model.
#
# Writes, per subject and contrast, a z map and an effect-size map. The effect
# maps are the ones to take to the group level.

# ---------------------------------------------------------------- config ----
BIDS_DIR = "/projects/bd2/bids"
FMRIPREP_DIR = "/projects/bd2/derivatives/fmriprep"
OUT_DIR = "/projects/bd2/derivatives/first_level_mid"
SPACE = "MNI152NLin2009cAsym"

SESSION = "ses-1"
TASK = "mid"
SMOOTHING_FWHM = 5.0          # mm; data are 2 mm isotropic
HRF_MODEL = "spm"
NOISE_MODEL = "ar1"

# fMRIPrep's cosine regressors ARE the high-pass filter, so drift_model is None
# below. Adding a second drift model on top would double-filter.
MOTION = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]
TISSUE = ["csf", "white_matter"]
# -----------------------------------------------------------------------------

# Events use "reward.high" etc. The dot breaks nilearn's contrast expression
# parser, so trial types are sanitised to reward_high on load.
CONTRASTS = {
    # reward anticipation vs implicit baseline, collapsed over magnitude
    "antReward":
        "0.5*ant_reward_high + 0.5*ant_reward_low",
    # positive reward feedback (a Hit on a win trial = money won) vs baseline
    "fbRewardPos":
        "0.5*rew_reward_high_Hit + 0.5*rew_reward_low_Hit",
    # the same two, each minus its neutral ($0) counterpart
    "antRewardMinusNeut":
        "0.5*ant_reward_high + 0.5*ant_reward_low - ant_reward_neut",
    "fbRewardPosMinusNeut":
        "0.5*rew_reward_high_Hit + 0.5*rew_reward_low_Hit - rew_reward_neut_Hit",
}


def sanitise(events):
    events = events.copy()
    events["trial_type"] = events["trial_type"].str.replace(".", "_", regex=False)
    return events


def run_inputs(sub, run, args):
    """Paths for one run, or None with a reason if anything is missing."""
    stem = "%s_%s_task-%s_run-%s" % (sub, SESSION, TASK, run)
    func = os.path.join(args.fmriprep, sub, SESSION, "func")
    bold = os.path.join(func, "%s_space-%s_desc-preproc_bold.nii.gz" % (stem, args.space))
    mask = os.path.join(func, "%s_space-%s_desc-brain_mask.nii.gz" % (stem, args.space))
    conf = os.path.join(func, "%s_desc-confounds_timeseries.tsv" % stem)
    events = os.path.join(args.bids, sub, SESSION, "func", "%s_events.tsv" % stem)
    sidecar = os.path.join(args.bids, sub, SESSION, "func", "%s_bold.json" % stem)
    for label, path in (("preproc bold", bold), ("brain mask", mask),
                        ("confounds", conf), ("events", events), ("sidecar", sidecar)):
        if not os.path.exists(path):
            return None, "missing %s (%s)" % (label, os.path.basename(path))
    return dict(bold=bold, mask=mask, conf=conf, events=events, sidecar=sidecar), None


def load_confounds(path):
    """Motion + tissue + fMRIPrep's cosine drift terms."""
    df = pd.read_csv(path, sep="\t")
    cols = [c for c in MOTION + TISSUE if c in df.columns]
    cols += [c for c in df.columns if c.startswith("cosine")]
    out = df[cols].copy()
    # First-row NaNs come from derivative-style regressors; a zero there is the
    # standard handling and keeps the design full rank.
    return out.fillna(0.0)


def fit_subject(sub, args):
    runs, skipped = [], []
    for run in ("01", "02"):
        info, why = run_inputs(sub, run, args)
        (skipped.append("run-%s: %s" % (run, why)) if info is None else runs.append(info))
    if not runs:
        return "SKIP", "; ".join(skipped)

    t_rs = {json.load(open(r["sidecar"]))["RepetitionTime"] for r in runs}
    if len(t_rs) != 1:
        return "SKIP", "runs disagree on RepetitionTime: %s" % sorted(t_rs)
    t_r = t_rs.pop()

    events = [sanitise(pd.read_csv(r["events"], sep="\t")) for r in runs]
    present = set().union(*(set(e["trial_type"]) for e in events))
    missing = [c for c, expr in CONTRASTS.items()
               for term in _terms(expr) if term not in present]
    if missing:
        return "SKIP", "conditions absent from events: %s" % sorted(set(missing))

    if args.dry_run:
        return "OK", "%d run(s), TR=%s, would fit and write %d contrasts" % (
            len(runs), t_r, len(CONTRASTS))

    model = FirstLevelModel(
        t_r=t_r,
        hrf_model=HRF_MODEL,
        drift_model=None,               # cosine regressors supply the high-pass
        smoothing_fwhm=args.fwhm,
        mask_img=runs[0]["mask"],
        noise_model=NOISE_MODEL,
        standardize=False,
        signal_scaling=0,               # percent signal change relative to run mean
        minimize_memory=False,
    )
    model.fit(
        run_imgs=[r["bold"] for r in runs],
        events=events,
        confounds=[load_confounds(r["conf"]) for r in runs],
    )

    out = os.path.join(args.out, sub)
    os.makedirs(out, exist_ok=True)
    for name, expr in CONTRASTS.items():
        for kind, suffix in (("z_score", "z"), ("effect_size", "effect")):
            img = model.compute_contrast(expr, output_type=kind)
            img.to_filename(os.path.join(
                out, "%s_task-%s_contrast-%s_stat-%s_statmap.nii.gz"
                % (sub, TASK, name, suffix)))
    return "OK", "%d run(s), TR=%s, %d contrasts written" % (len(runs), t_r, len(CONTRASTS))


def _terms(expr):
    """Condition names referenced by a contrast expression."""
    for chunk in expr.replace("-", "+").split("+"):
        chunk = chunk.strip()
        if "*" in chunk:
            chunk = chunk.split("*", 1)[1].strip()
        if chunk:
            yield chunk


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("subjects", nargs="*", help="subject labels (default: all in fmriprep dir)")
    ap.add_argument("--bids", default=BIDS_DIR)
    ap.add_argument("--fmriprep", default=FMRIPREP_DIR)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--space", default=SPACE)
    ap.add_argument("--fwhm", type=float, default=SMOOTHING_FWHM)
    ap.add_argument("--dry-run", action="store_true", help="check inputs, fit nothing")
    args = ap.parse_args()

    found = sorted(os.path.basename(p) for p in glob.glob(args.fmriprep + "/sub-*")
                   if os.path.isdir(p))
    if not found:
        sys.exit("no sub-* directories under %s -- is --fmriprep right?" % args.fmriprep)
    subs = (["sub-" + s.replace("sub-", "") for s in args.subjects]
            if args.subjects else found)

    failures = 0
    for sub in subs:
        try:
            status, detail = fit_subject(sub, args)
        except Exception as exc:                       # keep the batch going
            status, detail = "FAIL", "%s: %s" % (type(exc).__name__, exc)
        if status != "OK":
            failures += 1
        print("%-7s %-11s %s" % (status, sub, detail), flush=True)

    print("\n%d subject(s), %d not fit" % (len(subs), failures))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
