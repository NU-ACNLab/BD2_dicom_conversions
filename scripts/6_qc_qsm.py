import os
import sys

import nibabel as nib
import numpy as np

from qsm_common import (BIDS_DIR, EXPECTED_ECHOES, discover_subjects,
                        echo_files, read_te_ms, session_dir)

# QC for MEGRE (QSM source) data.
#
# Usage: python 6_qc_qsm.py [SES] [SUB ...] [--tsv PATH]
#
#   python 6_qc_qsm.py ses-1                 # every subject with that session
#   python 6_qc_qsm.py ses-1 rb860 rb926     # just these
#   python 6_qc_qsm.py ses-1 --tsv qc.tsv    # also write a table
#
# Checks the things that are cheap, automatable, and silent when they break --
# conversion and metadata integrity. It does NOT replace looking at the images:
# motion ghosting, phase singularities from failed coil combination, and
# background-field streaking all need eyes. See README for the visual pass.
#
# Exits 1 if any subject FAILs, so it can gate a batch run.

# Phase must be in radians before any QSM pipeline sees it (see
# 2b_phase_to_radians.py). Slack above pi for interpolation overshoot.
RADIAN_CEILING = 3.25

# Non-uniform echo spacing is legal -- pipelines that read per-echo TEs handle
# it -- but it breaks tools that assume a constant delta, so it is worth
# flagging. Tolerance in ms on the spread of successive TE differences.
DTE_SPREAD_TOL = 0.30

# Below this z-extent the volume is a slab, not a whole head, and background
# field removal will erode usable slices off both ends.
WHOLE_HEAD_Z_MM = 160.0

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"

TSV_COLUMNS = [
    "subject", "session", "n_mag", "n_phase", "te_first_ms", "te_last_ms",
    "dte_min_ms", "dte_max_ms", "dims", "voxel_mm", "z_extent_mm",
    "phase_min", "phase_max", "phase_units", "t2star_ms", "status",
]


def check_inventory(mag, phase, row):
    row["n_mag"], row["n_phase"] = len(mag), len(phase)
    detail = "%d mag + %d phase" % (len(mag), len(phase))

    if not mag and not phase:
        return FAIL, "no MEGRE files found"
    missing = sorted(set(mag) ^ set(phase))
    if missing:
        return FAIL, detail + "; echoes without a mag/phase pair: %s" % missing
    if len(mag) != EXPECTED_ECHOES:
        return WARN, detail + " (expected %d echoes)" % EXPECTED_ECHOES
    return PASS, detail


def check_echo_times(mag, phase, row):
    """TEs present, monotonic, and consistent between magnitude and phase."""
    echoes = sorted(mag)
    tes = [read_te_ms(mag[e]) for e in echoes]
    if any(t is None for t in tes):
        return FAIL, "missing EchoTime in one or more sidecars"

    row["te_first_ms"] = "%.3f" % tes[0]
    row["te_last_ms"] = "%.3f" % tes[-1]

    if any(b <= a for a, b in zip(tes, tes[1:])):
        return FAIL, "TEs not monotonic: %s" % [round(t, 3) for t in tes]

    for e in echoes:
        mag_te, phase_te = read_te_ms(mag[e]), read_te_ms(phase[e])
        if phase_te is None or abs(mag_te - phase_te) > 1e-6:
            return FAIL, "echo %d: mag TE %s != phase TE %s" % (e, mag_te, phase_te)

    span = "%.3f .. %.3f ms" % (tes[0], tes[-1])
    if len(tes) < 3:
        return PASS, span + ", monotonic"

    deltas = [b - a for a, b in zip(tes, tes[1:])]
    row["dte_min_ms"] = "%.3f" % min(deltas)
    row["dte_max_ms"] = "%.3f" % max(deltas)
    spread = max(deltas) - min(deltas)
    if spread > DTE_SPREAD_TOL:
        return WARN, (span + ", monotonic; NON-UNIFORM spacing "
                      "(dTE %.2f..%.2f ms, spread %.2f) -- confirm against the "
                      "protocol and use a pipeline that reads per-echo TEs"
                      % (min(deltas), max(deltas), spread))
    return PASS, span + ", monotonic, uniform dTE %.2f ms" % np.mean(deltas)


def check_geometry(mag, phase, row):
    """Shape, voxel size, and affine identical across every echo."""
    paths = [mag[e] for e in sorted(mag)] + [phase[e] for e in sorted(phase)]
    reference = nib.load(paths[0])
    shape, affine = reference.shape, reference.affine
    zooms = reference.header.get_zooms()[:3]

    row["dims"] = "x".join(str(s) for s in shape)
    row["voxel_mm"] = "x".join("%.3f" % z for z in zooms)
    z_extent = shape[2] * zooms[2]
    row["z_extent_mm"] = "%.1f" % z_extent

    for path in paths[1:]:
        img = nib.load(path)
        if img.shape != shape:
            return FAIL, "%s has shape %s, expected %s" % (
                os.path.basename(path), img.shape, shape)
        if not np.allclose(img.affine, affine, atol=1e-4):
            return FAIL, "%s affine differs from echo 1" % os.path.basename(path)

    detail = "%s @ %s mm, consistent across %d files" % (
        row["dims"], row["voxel_mm"], len(paths))
    if z_extent < WHOLE_HEAD_Z_MM:
        return WARN, (detail + "; z-extent %.0f mm is a slab -- background "
                      "field removal will erode both ends, confirm your ROIs "
                      "survive" % z_extent)
    return PASS, detail


def check_phase_units(phase, row):
    """The big one: phase must be radians, not raw scanner units."""
    if not phase:
        return FAIL, "no phase images"

    lo = hi = None
    for e in sorted(phase):
        data = nib.load(phase[e]).get_fdata()
        lo = float(data.min()) if lo is None else min(lo, float(data.min()))
        hi = float(data.max()) if hi is None else max(hi, float(data.max()))

    row["phase_min"], row["phase_max"] = "%.4f" % lo, "%.4f" % hi
    span = "range [%.3f, %.3f]" % (lo, hi)

    if -RADIAN_CEILING <= lo and hi <= RADIAN_CEILING:
        row["phase_units"] = "rad"
        # A near-zero span means the volume is flat, not that it is in radians.
        if hi - lo < 1.0:
            return WARN, span + " -- suspiciously narrow for wrapped phase"
        return PASS, span + " -- radians"

    row["phase_units"] = "scanner"
    return FAIL, (span + " -- NOT radians. QSM pipelines will not error, they "
                  "will produce a chi map wrong by ~%dx. Run "
                  "2b_phase_to_radians.py." % int(hi / np.pi))


def check_mag_decay(mag, row):
    """Signal must fall monotonically across the echo train.

    A rise means an echo is mis-sorted or duplicated -- the failure mode that
    a filename-based pipeline cannot otherwise see.
    """
    echoes = sorted(mag)
    first = nib.load(mag[echoes[0]]).get_fdata()

    # Crude tissue mask; adequate for a relative signal trend.
    mask = first > 0.2 * np.percentile(first, 99.5)
    if mask.sum() < 1000:
        return WARN, "could not build a tissue mask from echo 1"

    means = [float(nib.load(mag[e]).get_fdata()[mask].mean()) for e in echoes]

    tes = [read_te_ms(mag[e]) for e in echoes]
    if all(t is not None for t in tes) and means[-1] > 0:
        # Log-linear fit; slope is -1/T2*.
        slope = np.polyfit(tes, np.log(means), 1)[0]
        t2star = -1.0 / slope if slope < 0 else float("nan")
        row["t2star_ms"] = "%.1f" % t2star
        tail = ", T2* ~ %.0f ms" % t2star
    else:
        tail = ""

    rises = [echoes[i + 1] for i in range(len(means) - 1) if means[i + 1] > means[i]]
    detail = "%.1f -> %.1f%s" % (means[0], means[-1], tail)
    if rises:
        return FAIL, detail + "; signal RISES at echo(es) %s -- echoes likely " \
                              "mis-sorted or duplicated" % rises
    return PASS, "monotonic decay, " + detail


def qc_session(sub, ses):
    """Run every check for one subject/session. Returns (status, row)."""
    session = session_dir(sub, ses)
    row = dict.fromkeys(TSV_COLUMNS, "")
    row["subject"], row["session"] = "sub-" + sub, ses

    print("\nsub-%s / %s" % (sub, ses))

    mag, phase = echo_files(session, "mag"), echo_files(session, "phase")

    checks = [("inventory", lambda: check_inventory(mag, phase, row))]
    if mag and phase and not (set(mag) ^ set(phase)):
        checks += [
            ("echo times", lambda: check_echo_times(mag, phase, row)),
            ("geometry", lambda: check_geometry(mag, phase, row)),
            ("phase units", lambda: check_phase_units(phase, row)),
            ("mag decay", lambda: check_mag_decay(mag, row)),
        ]

    worst = PASS
    for label, check in checks:
        try:
            status, detail = check()
        except Exception as exc:                      # a broken file is a FAIL
            status, detail = FAIL, "%s: %s" % (type(exc).__name__, exc)
        print("  [%s] %-12s %s" % (status, label, detail))
        if status == FAIL or (status == WARN and worst == PASS):
            worst = status

    row["status"] = worst
    return worst, row


def write_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in TSV_COLUMNS) + "\n")
    print("\nwrote %s" % path)


def main():
    argv = sys.argv[1:]
    tsv_path = None
    if "--tsv" in argv:
        i = argv.index("--tsv")
        if i + 1 >= len(argv):
            sys.exit("--tsv needs a path")
        tsv_path = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]

    if not argv:
        sys.exit("usage: %s SES [SUB ...] [--tsv PATH]"
                 % os.path.basename(sys.argv[0]))

    ses, subs = argv[0], argv[1:]
    if not subs:
        subs = discover_subjects(ses)
        if not subs:
            sys.exit("no subjects with session %r under %s" % (ses, BIDS_DIR))

    results = [qc_session(sub, ses) for sub in subs]
    rows = [row for _, row in results]

    tally = {PASS: 0, WARN: 0, FAIL: 0}
    for status, _ in results:
        tally[status] += 1
    print("\n%d subject(s): %d pass, %d warn, %d fail"
          % (len(results), tally[PASS], tally[WARN], tally[FAIL]))

    if tsv_path:
        write_tsv(tsv_path, rows)

    if tally[FAIL]:
        sys.exit(1)


main()
