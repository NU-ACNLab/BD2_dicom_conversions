import glob
import json
import os
import re

# Shared helpers for the MEGRE/QSM scripts (6_qc_qsm.py, 7_qc_images.py).
#
# These two need to agree exactly on which files are MEGRE echoes and what TE
# each one carries -- if they ever disagree, the QC table and the QC images
# describe different data. Hence one copy, imported by both.

BIDS_DIR = os.environ.get(
    "BD2_ROOT", "/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged"
) + "/bids"

EXPECTED_ECHOES = 8


def session_dir(sub, ses):
    return os.path.join(BIDS_DIR, "sub-" + sub, ses)


def echo_files(session_path, part):
    """Map echo number -> .nii.gz path for 'mag' or 'phase'.

    Searches recursively, so it works both before 3_deface.sh sorts the
    session into anat/ and after.
    """
    pattern = os.path.join(session_path, "**", "*part-%s*MEGRE.nii.gz" % part)
    found = {}
    for path in glob.glob(pattern, recursive=True):
        match = re.search(r"_echo-(\d+)_", os.path.basename(path))
        if match:
            found[int(match.group(1))] = path
    return found


def read_te_ms(nii_path):
    """EchoTime from the sidecar, in ms. None if unreadable."""
    json_path = nii_path[: -len(".nii.gz")] + ".json"
    if not os.path.exists(json_path):
        return None
    with open(json_path) as fh:
        te = json.load(fh).get("EchoTime")
    return None if te is None else te * 1000.0


def discover_subjects(ses):
    """Subject labels (without the 'sub-' prefix) that have this session."""
    subs = []
    for path in sorted(glob.glob(os.path.join(BIDS_DIR, "sub-*"))):
        if os.path.isdir(os.path.join(path, ses)):
            subs.append(os.path.basename(path)[len("sub-"):])
    return subs
