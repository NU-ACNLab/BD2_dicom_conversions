import argparse
import csv
import os
import re
import sys

import pandas as pd

# PsychoPy behavioral output -> BIDS events.tsv, for MID and RL.
#
# Usage: python 8_behavioral_to_events.py [SUB ...] [options]
#
#   python 8_behavioral_to_events.py                      # everyone
#   python 8_behavioral_to_events.py rs1356 rb1217        # just these
#   python 8_behavioral_to_events.py --dry-run            # write nothing
#   python 8_behavioral_to_events.py --qc events_qc.csv   # QC table elsewhere
#
# Onsets in the PsychoPy csvs are already zeroed at each run's scanner trigger
# (MID run 1 and run 2 both put their first cue at ~5.0 s; every RL run starts
# its first cue at ~0.001 s), so nothing is re-referenced here.
#
# Nothing is guessed. A run is written only when exactly one file supplies a
# complete copy of it. Restarts, aborts, partial runs and unparseable folders
# all land in the QC table instead, with the candidate files named, for you to
# resolve by hand.
#
# Exits 1 if any subject/run was skipped, so it can gate a batch run.

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BEH_DIR = os.path.join(ROOT, "behavioral")
BIDS_DIR = os.path.join(ROOT, "bids")
# Subjects with behavioral data but no bids/sub-XX yet: their events land here
# so they can be dropped in after conversion, rather than being lost or
# creating bold-less func dirs that fail validation.
STAGING_DIR = os.path.join(BEH_DIR, "events_pending")
QC_PATH = os.path.join(BEH_DIR, "events_qc.csv")

SESSION = "ses-1"

# behavioral/<site dir>/<folder> -> subject label. Berkeley exports the whole
# study name as the folder; UCL and Swinburne already use the label.
SITE_DIRS = {
    "UC Berkeley BD2 Behavioral Data": ("UCB", re.compile(r"^Johnson_Berk_RandR_(\d+)$"), "rb"),
    "UCL BD2 Behavioral Data": ("UCL", re.compile(r"^rl(\d+)$"), "rl"),
    "Swinburne BD2 Behavioral Data": ("SUT", re.compile(r"^rs(\d+)$"), "rs"),
}

# MID: run 0 is the MRT target-duration calibration block -- it has no cue and
# no BOLD run, so only runs 1 and 2 become events files.
MID_RUNS = {1: "01", 2: "02"}
MID_TRIALS_PER_RUN = 36
# RL: three scanned runs, zero-indexed in the csv.
RL_RUNS = {0: "01", 1: "02", 2: "03"}
RL_TRIALS_PER_RUN = 50

MID_COLUMNS = [
    "run", "trial.type", "Cue.OnsetTime", "Tgt.OnsetTime", "trial.rt",
    "Fix_after_target.OnsetTime", "Fb.OnsetTime", "Fix_ITI.OnsetTime",
]
RL_COLUMNS = [
    "runs.thisN", "cueOnTime", "cueOffTime", "cueRespOnTime", "cueRespOffTime",
    "outcomeOnTime", "outcomeOffTime", "cue_resp.keys", "good_side",
    "outcome_image",
]

SIDES = ("left", "right")

QC_COLUMNS = ["subject", "site", "task", "run", "status", "n_trials", "detail", "source"]

OK, PARTIAL, SKIP = "ok", "PARTIAL", "SKIP"


def discover_subjects():
    """Map every parseable behavioral folder to its subject label.

    Returns (subjects, unparsed) where subjects is {label: (site, path)} and
    unparsed is a list of (site, folder) that did not match a site's ID pattern
    -- pilots and test IDs, which are reported and not processed.
    """
    subjects, unparsed = {}, []
    for site_dir, (site, pattern, prefix) in sorted(SITE_DIRS.items()):
        path = os.path.join(BEH_DIR, site_dir)
        if not os.path.isdir(path):
            continue
        for folder in sorted(os.listdir(path)):
            full = os.path.join(path, folder)
            if not os.path.isdir(full):
                continue
            match = pattern.match(folder)
            if match is None:
                unparsed.append((site, folder))
                continue
            subjects[prefix + match.group(1)] = (site, full)
    return subjects, unparsed


def find_task_files(subject_dir, task):
    """Every in-scanner csv for a task, anywhere under the subject folder.

    Sites nest differently -- Berkeley uses <folder>_MID/, UCL uses MID/,
    Swinburne keeps everything flat -- so this walks rather than globbing a
    fixed depth. Practice, out-of-scanner, backup and the per-run
    target_durs sidecars are all excluded.
    """
    hits = []
    for dirpath, dirnames, filenames in os.walk(subject_dir):
        dirnames[:] = [d for d in dirnames if d.lower() not in ("practice", "out-of-scanner")]
        for name in filenames:
            if not name.endswith(".csv") or "backup" in name:
                continue
            if task == "mid":
                if not name.startswith("MID") or "_fmri_" not in name:
                    continue
                if "target_durs" in name:
                    continue
            else:
                if "_RL_in-scanner_" not in name:
                    continue
            hits.append(os.path.join(dirpath, name))
    return sorted(hits)


def read_task_file(path, task):
    """Load one csv, or None if it is not a usable copy of this task.

    PsychoPy leaves behind header-only and aborted files with the same name
    pattern as good ones; a file missing the timing columns is one of those.
    """
    try:
        df = pd.read_csv(path, low_memory=False)
    except Exception:
        return None
    needed = MID_COLUMNS if task == "mid" else RL_COLUMNS
    if any(col not in df.columns for col in needed):
        return None
    run_col, onset_col = (("run", "Cue.OnsetTime") if task == "mid"
                          else ("runs.thisN", "cueOnTime"))
    df = df[df[onset_col].notna() & df[run_col].notna()]
    return df if len(df) else None


def candidates_by_run(files, task):
    """{run index: [(path, trial dataframe), ...]} across every file found."""
    run_col = "run" if task == "mid" else "runs.thisN"
    runs = RL_RUNS if task == "rl" else MID_RUNS
    found = {run: [] for run in runs}
    for path in files:
        df = read_task_file(path, task)
        if df is None:
            continue
        for run in runs:
            trials = df[df[run_col].astype(float) == float(run)]
            if len(trials):
                found[run].append((path, trials))
    return found


def emit(rows, dropped, onset, duration, label):
    """Append one event, unless its timing is missing.

    Recovered and crash-truncated files can lose the tail of a trial -- a
    response logged but the feedback onset never written. Those events cannot
    be placed, so they are dropped and counted rather than written as NaN.
    """
    if pd.isna(onset) or pd.isna(duration):
        dropped.append(label)
        return
    rows.append((float(onset), float(duration), label))


def mid_events(trials):
    """MID phase rows: anticipation, target window, button press, feedback.

    Durations are measured from the logged onsets rather than the nominal
    values, because the target duration is staircased per trial and the
    anticipation window jitters by ~0.5 s.
    """
    rows, dropped = [], []
    for _, t in trials.iterrows():
        cond = t["trial.type"]
        cue, tgt = t["Cue.OnsetTime"], t["Tgt.OnsetTime"]
        fix_after, fbk = t["Fix_after_target.OnsetTime"], t["Fb.OnsetTime"]
        iti, rt = t["Fix_ITI.OnsetTime"], t["trial.rt"]
        # A response of any kind is a hit -- verified against Tgt.ACCfeedback,
        # which reads Miss for every trial with no RT and Hit for every trial
        # with one. Tgt.ACC itself splits misses across 0 and 3.
        hit = pd.notna(rt)

        emit(rows, dropped, cue, tgt - cue, "ant_" + cond)
        emit(rows, dropped, tgt, fix_after - tgt, "motor_period")
        if hit:
            emit(rows, dropped, tgt + rt, 0.0, "motor_response")
        emit(rows, dropped, fbk, iti - fbk, "rew_%s_%s" % (cond, "Hit" if hit else "Miss"))
    return rows, dropped


def rl_key_sides(trials):
    """Recover which response key means left and which means right.

    Response boxes differ by site and changed within UCL mid-study -- 1/2, 0/1,
    6/7, 2/3 and b/y all appear in these files -- so the mapping is read out of
    each run rather than hardcoded. Reward only ever follows a choice of
    good_side, so among rewarded trials a key appears with exactly one side.

    Returns {key: side}, or None if the run does not pin the mapping down.
    """
    keys = set(trials["cue_resp.keys"].dropna().astype(str))
    if len(keys) != 2:
        return None

    rewarded = trials[trials["outcome_image"].astype(str).str.contains("coin", na=False)]
    mapping = {}
    for key, group in rewarded.groupby(rewarded["cue_resp.keys"].astype(str)):
        sides = set(group["good_side"].astype(str))
        if len(sides) != 1:
            return None  # a key won under both sides -- the invariant is broken
        mapping[key] = sides.pop()

    if len(mapping) == 2:
        return mapping if len(set(mapping.values())) == 2 else None
    if len(mapping) == 1:
        # One key never won. It must be the other side.
        known_key, known_side = next(iter(mapping.items()))
        if known_side not in SIDES:
            return None
        other_key = (keys - {known_key}).pop()
        mapping[other_key] = SIDES[1 - SIDES.index(known_side)]
        return mapping
    return None


def rl_events(trials, key_sides):
    """RL phase rows: decision, selection display, outcome.

    The inter-phase fixation is left as implicit baseline. Reward only ever
    follows a choice of good_side, so outcome_reward and outcome_noreward
    already carry the choice-correctness contrast.
    """
    rows, dropped = [], []
    for _, t in trials.iterrows():
        cue_on, cue_off = t["cueOnTime"], t["cueOffTime"]
        sel_on, sel_off = t["cueRespOnTime"], t["cueRespOffTime"]
        out_on, out_off = t["outcomeOnTime"], t["outcomeOffTime"]
        key = t["cue_resp.keys"]
        side = key_sides.get(str(key)) if pd.notna(key) else None
        image = str(t["outcome_image"])

        if side is None:
            decision, selection, outcome = "decision_noresp", "select_noresp", "outcome_noresp"
        else:
            optimal = "optimal" if side == str(t["good_side"]) else "nonoptimal"
            decision = "decision"
            selection = "select_" + optimal
            outcome = "outcome_reward" if "coin" in image else "outcome_noreward"

        emit(rows, dropped, cue_on, cue_off - cue_on, decision)
        emit(rows, dropped, sel_on, sel_off - sel_on, selection)
        emit(rows, dropped, out_on, out_off - out_on, outcome)
    return rows, dropped


def write_events(rows, path):
    df = pd.DataFrame(rows, columns=["onset", "duration", "trial_type"])
    df = df.sort_values("onset", kind="stable").reset_index(drop=True)
    df["onset"] = df["onset"].round(3)
    df["duration"] = df["duration"].round(3)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return len(df)


def destination(subject, task, run_label, staging):
    """bids/sub-XX/ses-1/func/ if that subject is converted, else staging."""
    name = "sub-%s_%s_task-%s_run-%s_events.tsv" % (subject, SESSION, task, run_label)
    func = os.path.join(BIDS_DIR, "sub-" + subject, SESSION, "func")
    if os.path.isdir(func):
        return os.path.join(func, name), "bids"
    return os.path.join(staging, name), "staging"


def process(subject, site, subject_dir, task, args, qc):
    """Build every complete run of one task for one subject."""
    task_label = "mid" if task == "mid" else "RL"
    expected = MID_TRIALS_PER_RUN if task == "mid" else RL_TRIALS_PER_RUN
    runs = MID_RUNS if task == "mid" else RL_RUNS

    files = find_task_files(subject_dir, task)
    if not files:
        qc.append([subject, site, task_label, "", SKIP, "", "no in-scanner file found", ""])
        return False

    found = candidates_by_run(files, task)
    clean = True
    for run, run_label in sorted(runs.items()):
        cands = found[run]
        complete = [(p, t) for p, t in cands if len(t) == expected]

        if not cands:
            qc.append([subject, site, task_label, run_label, SKIP, 0,
                       "run absent from every file", rel(files)])
            clean = False
            continue
        if not complete:
            counts = ", ".join("%s=%d" % (os.path.basename(p), len(t)) for p, t in cands)
            qc.append([subject, site, task_label, run_label, SKIP,
                       max(len(t) for _, t in cands),
                       "partial run, expected %d trials (%s)" % (expected, counts),
                       rel([p for p, _ in cands])])
            clean = False
            continue
        if len(complete) > 1:
            qc.append([subject, site, task_label, run_label, SKIP, expected,
                       "%d files each hold a complete run -- restart, resolve by hand"
                       % len(complete),
                       rel([p for p, _ in complete])])
            clean = False
            continue

        source, trials = complete[0]
        if task == "mid":
            rows, dropped = mid_events(trials)
        else:
            key_sides = rl_key_sides(trials)
            if key_sides is None:
                qc.append([subject, site, task_label, run_label, SKIP, len(trials),
                           "cannot tell which response key is left vs right",
                           rel([source])])
                clean = False
                continue
            rows, dropped = rl_events(trials, key_sides)

        path, where = destination(subject, task_label, run_label, args.staging)
        if not args.dry_run:
            write_events(rows, path)

        detail = "%d events -> %s" % (len(rows), where)
        status = OK
        if dropped:
            status = PARTIAL
            detail += "; %d event(s) dropped, timing never logged (%s)" % (
                len(dropped), ", ".join(sorted(set(dropped))))
            clean = False
        qc.append([subject, site, task_label, run_label, status, len(trials), detail,
                   rel([source])])
    return clean


def rel(paths):
    return "; ".join(os.path.relpath(p, ROOT) for p in paths)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("subjects", nargs="*", help="subject labels, e.g. rs1356 (default: all)")
    ap.add_argument("--task", choices=["mid", "rl", "both"], default="both")
    ap.add_argument("--staging", default=STAGING_DIR,
                    help="where events for subjects without a bids/ dir go")
    ap.add_argument("--qc", default=QC_PATH, help="path for the QC table")
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    args = ap.parse_args()

    subjects, unparsed = discover_subjects()
    if args.subjects:
        missing = [s for s in args.subjects if s not in subjects]
        if missing:
            sys.exit("no behavioral folder for: " + ", ".join(missing))
        subjects = {s: subjects[s] for s in args.subjects}

    qc = []
    for site, folder in unparsed:
        qc.append(["", site, "", "", SKIP, "", "folder is not a study ID (pilot/test)", folder])

    tasks = ["mid", "rl"] if args.task == "both" else [args.task]
    clean = True
    for subject, (site, subject_dir) in sorted(subjects.items()):
        for task in tasks:
            clean &= process(subject, site, subject_dir, task, args, qc)

    if not args.dry_run:
        os.makedirs(os.path.dirname(args.qc), exist_ok=True)
        with open(args.qc, "w", newline="") as fo:
            writer = csv.writer(fo)
            writer.writerow(QC_COLUMNS)
            writer.writerows(qc)

    wrote = [r for r in qc if r[4] in (OK, PARTIAL)]
    flagged = [r for r in qc if r[4] != OK]
    to_bids = sum(1 for r in wrote if "-> bids" in r[6])
    print("%d run(s) written (%d into bids/, %d staged), %d skipped%s"
          % (len(wrote), to_bids, len(wrote) - to_bids,
             sum(1 for r in qc if r[4] == SKIP),
             " -- DRY RUN, nothing written" if args.dry_run else ""))
    for row in flagged:
        print("  %-7s %-8s %-4s run-%-2s  %s"
              % (row[4], row[0] or row[7], row[2], row[3], row[6]))
    if not args.dry_run:
        print("QC table: " + os.path.relpath(args.qc, ROOT))
    sys.exit(0 if clean else 1)


if __name__ == "__main__":
    main()
