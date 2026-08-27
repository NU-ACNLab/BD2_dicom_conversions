import os
import re
import sys
import glob

import pydicom

# Identify the site and participant of every DICOM folder from its headers.
#
# Usage: python 0_identify_dicoms.py [--tsv PATH] [--root DIR]
#
# Folder names under dicoms/uncompressed are inconsistent -- they arrive as
# whatever the exporting site happened to call them ("1-845", "P00<nnnn>",
# "rs1403_MR_1", "rs1403_MR_1 2"). The scanner headers are not inconsistent, so
# this reads the site and participant from the data itself and reports where
# the folder name disagrees.
#
# Read-only. It renames nothing and converts nothing; it prints a table and, if
# asked, writes a TSV. Use its output to decide what to rename by hand.

ROOT = os.environ.get(
    "BD2_ROOT", "/Users/katharinaseitz/Documents/BD2/BD2_dicom_conversions-merged"
)

# InstitutionName (0008,0080) is the primary key -- it is written by the
# scanner and none of the three sites share a value. DeviceSerialNumber is
# carried as an independent cross-check: if the two ever disagree, something is
# wrong with the export and the subject should not be converted blind.
SITE_BY_INSTITUTION = {
    "UC Berkeley": "UCB",
    "University College": "UCL",
    "Swinburne University": "SUT",
}

SITE_BY_SERIAL = {
    "166319": "UCB",
    "166108": "UCL",
    "167089": "SUT",
}

# How each site encodes the study participant number, and what we call it.
#
# Labels are r + a site letter: rb Berkeley, rl London, rs Swinburne.
#
#   UCB  PatientID   "Johnson_Berk_RandR_<n>" -> rb860
#   SUT  PatientName "2514xRNR_<n>"          -> rs1253
#        (SUT's PatientID is a scanner-local number like <6-digit local id> -- NOT the
#         study ID. Reading PatientID here would silently produce wrong labels.)
#
# UCL has NO header rule, deliberately. Its scans carry no study participant
# number anywhere: PatientID is UCL's own sequential scan accession (P00<nnnn>,
# P00<nnnn>, ...), PatientName is a deidentified date-time stamp, and a search of
# every element of every series for the known study numbers finds nothing.
# Deriving a label from the P-number gives a plausible but wrong subject ID,
# which is worse than none -- so UCL falls back to FOLDER_LABEL_SITES below.


def _berkeley(number):
    return "rb" + number


def _swinburne(number):
    return "rs" + number


PARTICIPANT_RULES = {
    "UCB": ("PatientID", r"Johnson_Berk_RandR_(\d+)$", _berkeley),
    "SUT": ("PatientName", r"2514xRNR_(\d+)$", _swinburne),
}

# Sites whose study number exists only in the folder name, and the label prefix
# to give it. The folder name is the operator's record of which participant a
# scan belongs to; rename the folder and the label follows.
FOLDER_LABEL_SITES = {"UCL": "rl"}

# Label prefix -> site. Lets a subject already converted into bids/ be assigned
# to its site without re-deriving it from the DICOMs, which matters when the
# source folder was never renamed.
LABEL_PREFIX_SITE = {"rb": "UCB", "rl": "UCL", "rs": "SUT"}

# Folder spellings accepted for those sites, each yielding the study number:
#   rl845 / RL845   already labelled
#   1-845           the exported "<batch>-<number>" form
#   845             bare number
# Anything else -- notably UCL's own P00<nnnn> -- carries no study number and is
# reported so the folder can be renamed.
FOLDER_LABEL_PATTERNS = (
    r"^(?:%(p)s)0*(\d+)$",
    r"^\d+-0*(\d+)$",
    r"^0*(\d+)$",
)


def label_from_folder(folder, prefix, flags):
    for pattern in FOLDER_LABEL_PATTERNS:
        match = re.match(pattern % {"p": prefix}, folder, re.IGNORECASE)
        if match:
            return "%s%d" % (prefix, int(match.group(1)))
    flags.append("NEEDS RENAME: no study number in folder name %r -- rename it "
                 "to %s<number> (the headers do not carry one)" % (folder, prefix))
    return None

TSV_COLUMNS = ["folder", "site", "participant", "study_date", "n_series",
               "n_files", "patient_id", "patient_name", "scanner", "software",
               "study_uid", "flags"]


def first_dicom(folder):
    """Path of the first readable file under folder, or None."""
    for path in sorted(glob.iglob(os.path.join(glob.escape(folder), "**", "*"),
                                  recursive=True)):
        if os.path.isfile(path) and not os.path.basename(path).startswith("."):
            return path
    return None


def read_header(path):
    try:
        return pydicom.dcmread(path, stop_before_pixels=True, force=True)
    except Exception:
        return None


def tag(ds, name):
    value = getattr(ds, name, None)
    return "" if value is None else str(value).strip()


def identify_site(ds, flags):
    """Site from InstitutionName, cross-checked against DeviceSerialNumber."""
    institution = tag(ds, "InstitutionName")
    serial = tag(ds, "DeviceSerialNumber")

    by_name = SITE_BY_INSTITUTION.get(institution)
    by_serial = SITE_BY_SERIAL.get(serial)

    if by_name and by_serial and by_name != by_serial:
        flags.append("SITE CONFLICT: institution=%s serial=%s" % (by_name, by_serial))
        return None
    if by_name:
        return by_name
    if by_serial:
        flags.append("site from serial only (institution %r unknown)" % institution)
        return by_serial

    flags.append("UNKNOWN SITE: institution=%r serial=%r" % (institution, serial))
    return None


def identify_participant(ds, site, folder, flags):
    """Study participant label, from the headers where they carry it and from
    the folder name where they do not."""
    if site in FOLDER_LABEL_SITES:
        return label_from_folder(folder, FOLDER_LABEL_SITES[site], flags)

    tag_name, pattern, label_for = PARTICIPANT_RULES[site]
    raw = tag(ds, tag_name)

    match = re.search(pattern, raw)
    if not match:
        flags.append("participant unparsed: %s=%r does not match %s for %s"
                     % (tag_name, raw, pattern, site))
        return None
    return label_for(match.group(1))


def scan_folder(folder_path):
    """Identify one DICOM folder. Returns a row dict."""
    folder = os.path.basename(folder_path)
    flags = []
    row = dict.fromkeys(TSV_COLUMNS, "")
    row["folder"] = folder

    series = [d for d in sorted(os.listdir(folder_path))
              if os.path.isdir(os.path.join(folder_path, d))]
    files = [f for f in glob.iglob(os.path.join(glob.escape(folder_path), "**", "*"),
                                   recursive=True) if os.path.isfile(f)]
    row["n_series"] = len(series)
    row["n_files"] = len(files)

    path = first_dicom(folder_path)
    if path is None:
        row["flags"] = "EMPTY: no files"
        return row

    ds = read_header(path)
    if ds is None:
        row["flags"] = "UNREADABLE: not a DICOM?"
        return row

    row["patient_id"] = tag(ds, "PatientID")
    row["patient_name"] = tag(ds, "PatientName")
    row["study_date"] = tag(ds, "StudyDate")
    row["study_uid"] = tag(ds, "StudyInstanceUID")
    row["scanner"] = tag(ds, "ManufacturerModelName")
    row["software"] = tag(ds, "SoftwareVersions")

    site = identify_site(ds, flags)
    row["site"] = site or "?"

    participant = identify_participant(ds, site, folder, flags) if site else None
    row["participant"] = participant or "?"

    # Report a folder whose name differs from the label rather than acting on
    # it -- renaming is the operator's call, and the pipeline passes the source
    # folder through with --src, so a mismatch is not itself a problem.
    if participant and folder != participant:
        flags.append("folder %r -> label %s" % (folder, participant))

    row["flags"] = "; ".join(flags)
    return row


def flag_duplicates(rows):
    """Mark folders that share a StudyInstanceUID -- the same scan, twice."""
    by_uid = {}
    for row in rows:
        if row["study_uid"]:
            by_uid.setdefault(row["study_uid"], []).append(row)

    for uid, group in by_uid.items():
        if len(group) < 2:
            continue
        names = [r["folder"] for r in group]
        for row in group:
            others = [n for n in names if n != row["folder"]]
            note = "DUPLICATE STUDY: same scan as %s" % ", ".join(others)
            row["flags"] = (row["flags"] + "; " + note) if row["flags"] else note


def read_existing_participants(path):
    """Existing participants.tsv as {participant_id: {column: value}}.

    source_id, scanner, software_version and study_date exist only in the DICOM
    headers. Once a converted subject's DICOMs are deleted -- a reasonable thing
    to do, they are bulky and the NIfTIs are the deliverable -- nothing else
    records them, so they are carried forward from the file we already wrote
    rather than being rewritten as n/a.
    """
    if not os.path.exists(path):
        return {}
    previous = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(header):
                continue          # malformed row: ignore rather than guess
            entry = dict(zip(header, values))
            if entry.get("participant_id"):
                previous[entry["participant_id"]] = entry
    return previous


def write_participants(root, rows):
    """Write bids/participants.tsv for subjects that have been converted.

    Only subjects with a bids/sub-LABEL directory are listed -- naming a
    subject that is not in the dataset is a validator error, and the DICOM
    roster runs ahead of what has actually been converted.

    site is the column that makes the single-dataset layout work: it is the
    only place the three cohorts are distinguishable once the files are
    renamed. source_id keeps each label traceable to the site's own records.
    """
    bids_dir = os.path.join(root, "bids")
    path = os.path.join(bids_dir, "participants.tsv")
    columns = ["participant_id", "site", "source_id", "scanner",
               "software_version", "study_date"]
    previous = read_existing_participants(path)

    # Driven by what is in bids/, not by the DICOM manifest. A subject
    # converted with an explicit label (4_single_sub.sh UCL rl169 --src ...)
    # has no identifiable manifest row while its source folder keeps the site's
    # own name, and omitting a subject that exists in the dataset is a
    # validator error. The label prefix always resolves the site.
    by_label = {}
    for row in rows:
        if row["participant"] != "?":
            by_label.setdefault(row["participant"], row)

    listed, carried = [], []
    for entry in sorted(glob.glob(os.path.join(bids_dir, "sub-*"))):
        if not os.path.isdir(entry):
            continue
        label = os.path.basename(entry)[len("sub-"):]
        row = by_label.get(label)

        participant_id = "sub-" + label
        was = previous.get(participant_id, {})

        def value(column, fresh):
            """DICOM header first, then what we recorded last time, then n/a."""
            if fresh:
                return fresh
            kept = was.get(column, "")
            if kept and kept != "n/a":
                carried.append("%s/%s" % (participant_id, column))
                return kept
            return "n/a"

        site = (row["site"] if row else "") or LABEL_PREFIX_SITE.get(label[:2], "")
        site = value("site", site)
        if site == "n/a":
            print("  sub-%s: unrecognised label prefix, site left as n/a" % label)

        listed.append({
            "participant_id": participant_id,
            "site": site,
            "source_id": value("source_id", row["patient_id"] if row else ""),
            "scanner": value("scanner", row["scanner"] if row else ""),
            "software_version": value("software_version", row["software"] if row else ""),
            "study_date": value("study_date", row["study_date"] if row else ""),
        })

    if not listed:
        print("\nno converted subjects yet -- participants.tsv not written")
        return

    if carried:
        subs = sorted({c.split("/")[0] for c in carried})
        print("\n  kept %d value(s) from the existing participants.tsv for %d "
              "subject(s) whose DICOMs are no longer present:\n    %s"
              % (len(carried), len(subs), ", ".join(subs)))

    listed.sort(key=lambda r: r["participant_id"])
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for entry in listed:
            fh.write("\t".join(entry[c] for c in columns) + "\n")

    sidecar = os.path.join(bids_dir, "participants.json")
    if not os.path.exists(sidecar):
        import json
        with open(sidecar, "w") as fh:
            json.dump({
                "site": {
                    "Description": "Acquisition site",
                    "Levels": {
                        "UCB": "University of California, Berkeley",
                        "UCL": "University College London",
                        "SUT": "Swinburne University of Technology",
                    },
                },
                "source_id": {
                    "Description": "Participant identifier in the acquiring "
                                   "site's own records (DICOM PatientID)",
                },
                "scanner": {"Description": "Scanner model (DICOM ManufacturerModelName)"},
                "software_version": {"Description": "Scanner software (DICOM SoftwareVersions)"},
                "study_date": {"Description": "Acquisition date, YYYYMMDD"},
            }, fh, indent=4)
            fh.write("\n")

    print("\nwrote %s (%d subject(s))" % (path, len(listed)))


def main():
    argv = sys.argv[1:]
    root = ROOT
    tsv_path = None
    participants = "--participants" in argv
    argv = [a for a in argv if a != "--participants"]

    for flag, setter in (("--tsv", "tsv"), ("--root", "root")):
        if flag in argv:
            i = argv.index(flag)
            if i + 1 >= len(argv):
                sys.exit("%s needs a value" % flag)
            if setter == "tsv":
                tsv_path = argv[i + 1]
            else:
                root = argv[i + 1]
            argv = argv[:i] + argv[i + 2:]

    if argv:
        sys.exit("usage: %s [--tsv PATH] [--participants] [--root DIR]"
                 % os.path.basename(sys.argv[0]))

    dicom_root = os.path.join(root, "dicoms", "uncompressed")
    if not os.path.isdir(dicom_root):
        sys.exit("no such directory: %s" % dicom_root)

    folders = [os.path.join(dicom_root, d) for d in sorted(os.listdir(dicom_root))
               if os.path.isdir(os.path.join(dicom_root, d))]
    if not folders:
        sys.exit("no subject folders under %s" % dicom_root)

    rows = [scan_folder(f) for f in folders]
    flag_duplicates(rows)

    width = max(len(r["folder"]) for r in rows)
    print("%-*s  %-4s  %-10s  %-8s  %5s  %6s" %
          (width, "FOLDER", "SITE", "PARTICIPANT", "DATE", "SER", "FILES"))
    for row in rows:
        print("%-*s  %-4s  %-10s  %-8s  %5s  %6s" %
              (width, row["folder"], row["site"], row["participant"],
               row["study_date"], row["n_series"], row["n_files"]))
        if row["flags"]:
            for note in row["flags"].split("; "):
                print("%s  ^ %s" % (" " * width, note))

    by_site = {}
    for row in rows:
        by_site[row["site"]] = by_site.get(row["site"], 0) + 1
    print("\n%d folder(s): %s" % (
        len(rows), ", ".join("%s=%d" % kv for kv in sorted(by_site.items()))))

    unresolved = [r for r in rows if r["participant"] == "?" or r["site"] == "?"]
    if unresolved:
        print("%d folder(s) could not be identified" % len(unresolved))

    if tsv_path:
        with open(tsv_path, "w") as fh:
            fh.write("\t".join(TSV_COLUMNS) + "\n")
            for row in rows:
                fh.write("\t".join(str(row[c]) for c in TSV_COLUMNS) + "\n")
        print("wrote %s" % tsv_path)

    if participants:
        write_participants(root, rows)


main()
