import os
import glob
import sys

# Merged UCB / UCL / SUT renamer.
#
# Usage: python 2_nifti_to_bids_naming.py SITE SUB SES [--dry-run]
#
# All three sites run the same pipeline; they differ only in what the scanner
# calls each series. Those differences live in SITES below -- nothing else in
# this file is site-specific.

BIDS_DIR = os.environ.get(
    "BD2_ROOT", "/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged"
) + "/bids"

# dcm2niix writes .nii.gz alongside sidecars; anything not listed here falls
# back to "everything after the first dot".
KNOWN_EXTS = (".nii.gz", ".nii", ".json", ".bval", ".bvec")


def get_ext(path):
    """Return the extension of a dcm2niix output, e.g. 'nii.gz' or 'json'.

    Replaces the old per-site file.split('.', 1)/split('.', 2) calls, which
    depended on how many periods a site's PatientName happened to contain and
    disagreed between sections of the same file.
    """
    base = os.path.basename(path)
    for ext in KNOWN_EXTS:
        if base.endswith(ext):
            return ext[1:]
    return base.split(".", 1)[1] if "." in base else ""


# Siemens writes derived maps alongside the diffusion series, named after it:
# "Axial MB DTI AP (MSV21)_ADC", "..._TRACEW", "..._FA", "..._ColFA". They are
# 3D scalar images, not diffusion data. "*DTI_AP*" matches all of them, and
# because a "one" rule renames every match to the same target they overwrite
# each other -- sub-rl845 ended up with the TRACEW map sitting in dwi/ as
# dir-AP_dwi, a single 3D volume where the real series has 7. UCL exports these
# maps; UCB and SUT do not, which is why only UCL was hit.
DERIVED_SUFFIXES = ("_ADC", "_TRACEW", "_FA", "_ColFA")


def is_derived(name):
    """True for a scanner-derived map that must not be renamed as raw data."""
    stem = os.path.basename(name).split("--")[0]
    return stem.endswith(DERIVED_SUFFIXES)


def _dwi_rules():
    # dir before run, per the BIDS dwi filename template -- "run-1_dir-AP_dwi"
    # is rejected by the validator. run is dropped rather than reordered:
    # there is one DTI pair, and the AP series is a b0-only reverse-phase scan
    # for topup, not a second run.
    return [
        ("one", "*DTI_AP*", "dir-AP_dwi"),
        ("one", "*DTI_PA*", "dir-PA_dwi"),
    ]


def _qsm_rules():
    """Tissue iron mapping: 8 magnitude echoes then 8 phase echoes.

    Named with the BIDS MEGRE suffix (multi-echo gradient echo), which is what
    QSM source data is -- not task fMRI. The old task-echo_flip-16_*_bold names
    described these as BOLD runs of a task that doesn't exist.
    """
    rules = []
    for i in range(1, 9):
        rules.append(("one", "*QSM_8echo*e%d.*" % i,
                      "echo-%d_part-mag_MEGRE" % i))
    for i in range(1, 9):
        rules.append(("one", "*QSM_8echo*e%d_ph.*" % i,
                      "echo-%d_part-phase_MEGRE" % i))
    return rules


def _rl_run_rules(template):
    """RL bold runs. UCB/SUT use 'task-RLCAT_run-0N', UCL uses 'task_RLCAT_run-0N'."""
    return [("one", "*" + template % n + "*", "task-RL_run-0%d_bold" % n)
            for n in (1, 2, 3)]


# Rules are applied in list order, and order matters: a file renamed by an
# earlier rule is no longer visible to a later one. The order below reproduces
# each branch's original ordering.
#
# Rule forms:
#   ("one", glob_pattern, bids_suffix)  -> rename every match to that suffix
#   ("t1w", glob_pattern, None)         -> T1w; 1st match becomes run-1, 2nd run-2, ...

SITES = {
    "UCB": [
        # AP_REV is the reverse-direction fieldmap; BIDS calls that dir-PA.
        ("one", "*--SpinEchoFieldmap_MID_AP_REV--*", "acq-mid_dir-PA_epi"),
        ("one", "*--SpinEchoFieldmap_MID_AP--*", "acq-mid_dir-AP_epi"),
        ("one", "*MID_run-1--*", "task-mid_run-01_bold"),
        ("one", "*MID_run-2--*", "task-mid_run-02_bold"),
        ("t1w", "*t1_mpg_sag*", None),
        ("one", "*--SpinEchoFieldmap_RLCAT_AP_REV--*", "acq-RL_dir-PA_epi"),
        ("one", "*--SpinEchoFieldmap_RLCAT_AP--*", "acq-RL_dir-AP_epi"),
    ] + _rl_run_rules("task-RLCAT_run-0%d") + _dwi_rules() + _qsm_rules(),

    "UCL": [
        ("one", "*SpinEchoFieldmap_MID_REV*", "acq-mid_dir-PA_epi"),
        ("one", "*SpinEchoFieldmap_MID_AP*", "acq-mid_dir-AP_epi"),
        # FIXED: the UCL branch had these two the other way round -- REV to
        # dir-AP -- which contradicted UCB, SUT, and UCL's own MID fieldmaps.
        # The sidecars settle it: at every site the _AP series carries
        # PhaseEncodingDirection j- and the _REV series carries j, for both
        # tasks. sub-rl845 was the only subject affected.
        # The REV rule must stay first: "*RLCAT_AP*" is unanchored and would
        # otherwise also swallow any *_AP_REV files.
        ("one", "*RLCAT_REV*", "acq-RL_dir-PA_epi"),
        ("one", "*RLCAT_AP*", "acq-RL_dir-AP_epi"),
        ("one", "*MID_run-1*", "task-mid_run-01_bold"),
        ("one", "*MID_run-2*", "task-mid_run-02_bold"),
        ("t1w", "*t1_mprage_sag*", None),
    ] + _rl_run_rules("task_RLCAT_run-0%d") + _dwi_rules() + _qsm_rules(),

    "SUT": [
        ("one", "*fmap_cmrr_mbep2d_se_PA--*", "acq-mid_dir-PA_epi"),
        ("one", "*fmap_cmrr_mbep2d_se_AP--*", "acq-mid_dir-AP_epi"),
        ("one", "*fmap_SpinEcho_RLCAT_AP--*", "acq-RL_dir-AP_epi"),
        ("one", "*fmap_SpinEcho_RLCAT_AP_REV--*", "acq-RL_dir-PA_epi"),
        ("one", "*MID_1*", "task-mid_run-01_bold"),
        ("one", "*MID_2*", "task-mid_run-02_bold"),
        ("t1w", "*t1_mprage_sag_p2_iso--*", None),
    ] + _rl_run_rules("task-RLCAT_run-0%d") + _dwi_rules() + _qsm_rules(),
}


def make_bids_dirs(sub, ses, bids_dir):
    for sub_dir in ("anat", "func", "fmap", "dwi"):
        os.makedirs(os.path.join(bids_dir, "sub-" + sub, ses, sub_dir), exist_ok=True)


def do_rename(old, new, dry_run):
    # A "one" rule renames every match to the same target, so two matches
    # silently overwrite each other and the survivor looks like a clean
    # conversion. That is how a derived TRACEW map ended up in sub-rl845's
    # dwi/. Refuse instead, and say which two files collided.
    if os.path.exists(new) and os.path.abspath(old) != os.path.abspath(new):
        raise SystemExit(
            "COLLISION: %s and the already-renamed %s both map to %s.\n"
            "Two source series match the same rule. Narrow the glob in SITES "
            "(or add the series to DERIVED_SUFFIXES) rather than letting one "
            "overwrite the other." % (os.path.basename(old),
                                      os.path.basename(new),
                                      os.path.basename(new)))
    print("%s -> %s" % (old, os.path.basename(new)))
    if not dry_run:
        os.rename(old, new)


def matches(directory, pattern):
    """Files matching a rule's glob, with scanner-derived maps filtered out."""
    return [f for f in sorted(glob.glob(directory + pattern), key=str)
            if not is_derived(f)]


def apply_rule(rule, sub, ses, directory, dry_run):
    kind, pattern, suffix = rule

    if kind == "one":
        for file in matches(directory, pattern):
            new_name = os.path.join(
                directory, "sub-%s_%s_%s.%s" % (sub, ses, suffix, get_ext(file)))
            do_rename(file, new_name, dry_run)
        return

    if kind == "t1w":
        # .json and .nii.gz are numbered independently so a subject's sidecar
        # and image get the same run number. sorted() keeps that assignment
        # stable -- the originals used raw glob order, which is arbitrary.
        for ext_glob in (".json", ".nii.gz"):
            files = matches(directory, pattern + ext_glob)
            for run, file in enumerate(files, start=1):
                new_name = os.path.join(
                    directory,
                    "sub-%s_%s_run-%d_T1w.%s" % (sub, ses, run, get_ext(file)))
                do_rename(file, new_name, dry_run)
        return

    raise ValueError("unknown rule kind: %s" % kind)


def rename_partic(site, sub, ses, bids_dir, dry_run):
    directory = bids_dir + "/sub-" + sub + "/" + ses + "/"
    print("renaming %s (site %s) in %s" % (sub, site, directory))
    for rule in SITES[site]:
        apply_rule(rule, sub, ses, directory, dry_run)


def main():
    args = [a for a in sys.argv[1:] if a != "--dry-run"]
    dry_run = "--dry-run" in sys.argv

    if len(args) != 3:
        sys.exit("usage: %s SITE SUB SES [--dry-run]   (SITE: %s)"
                 % (os.path.basename(sys.argv[0]), "|".join(sorted(SITES))))

    site, sub, ses = args
    site = site.upper()
    if site not in SITES:
        sys.exit("unknown site %r -- expected one of: %s"
                 % (site, ", ".join(sorted(SITES))))

    make_bids_dirs(sub, ses, BIDS_DIR)
    rename_partic(site, sub, ses, BIDS_DIR, dry_run)


if __name__ == "__main__":
    main()
