import json
import os
import sys

import matplotlib
matplotlib.use("Agg")            # must precede pyplot; there is no display here
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib            # noqa: E402
import numpy as np               # noqa: E402

from qsm_common import (BIDS_DIR, discover_subjects, echo_files,  # noqa: E402
                        read_te_ms, session_dir)

# Render QC images for MEGRE (QSM source) data.
#
# Usage: python 7_qc_images.py SES [SUB ...] [--outdir DIR]
#
#   python 7_qc_images.py ses-1                  # every subject with that session
#   python 7_qc_images.py ses-1 rb860            # just one
#
# Writes into derivatives/qc/sub-X/SES/ :
#
#   *_mag-echo1.png     anatomy, coverage, motion ghosting
#   *_mag-echoN.png     SNR and dropout at the longest TE
#   *_phase-echoN.png   THE one to look at -- see below
#   *_coverage.png      sagittal/coronal, slab position vs your ROIs
#
# ...plus a cross-subject contact sheet at the top of derivatives/qc/ so a
# batch can be reviewed in one glance rather than subject by subject.
#
# The phase montage uses the LAST echo deliberately. Phase accrues with TE, so
# a broken echo train still looks fine at 3 ms. What you are hunting for is an
# OPEN-ENDED FRINGE LINE: a black/white wrap boundary that starts inside the
# brain and stops mid-parenchyma instead of running out to the brain edge.
# That is a phase singularity from failed coil combination. Every unwrapping
# algorithm fails around it, and the resulting chi map looks plausible with a
# garbage region inside it. No automated metric in 6_qc_qsm.py catches it.

QC_DIR = os.path.join(os.path.dirname(BIDS_DIR), "derivatives", "qc")

N_SLICES = 12                 # 3x4 montage
GRID = (3, 4)
DPI = 95

# Fraction of the magnitude 99.5th percentile above which a voxel counts as
# tissue. Used to find the brain's z-extent so slices are sampled across the
# head rather than across the FOV (which wastes panels on empty air).
TISSUE_FRAC = 0.20


def brain_z_range(vol):
    """(z_lo, z_hi) covering the slices that actually contain signal."""
    thr = TISSUE_FRAC * np.percentile(vol, 99.5)
    per_slice = (vol > thr).sum(axis=(0, 1))
    live = np.nonzero(per_slice > 0.002 * vol.shape[0] * vol.shape[1])[0]
    if live.size < 2:
        return 0, vol.shape[2] - 1
    return int(live[0]), int(live[-1])


def display_limits(vol, phase):
    """Window for imshow. Phase is symmetric; magnitude clips the bright tail.

    Symmetric limits mean the phase panel renders identically whether the data
    is still in scanner units or has been through 2b_phase_to_radians.py.
    """
    if phase:
        lim = float(np.abs(vol).max()) or 1.0
        return dict(cmap="gray", vmin=-lim, vmax=lim)
    return dict(cmap="gray", vmin=0, vmax=float(np.percentile(vol, 99.5)))


def montage(vol, title, out_path, phase=False):
    """Axial montage sampled across the brain's z-extent."""
    lo, hi = brain_z_range(vol)
    zs = np.linspace(lo, hi, N_SLICES + 2)[1:-1]      # trim the end slices
    kw = display_limits(vol, phase)

    fig, axes = plt.subplots(*GRID, figsize=(13, 10.5), facecolor="k")
    for ax, z in zip(axes.ravel(), zs):
        ax.imshow(np.rot90(vol[:, :, int(round(z))]), **kw)
        # Boxed, or the label vanishes against a bright corner in phase panels.
        ax.text(3, 16, "z=%d" % int(round(z)), color="yellow", fontsize=8,
                va="top", bbox=dict(facecolor="black", edgecolor="none",
                                    alpha=0.65, pad=1.5))
        ax.set_axis_off()
    for ax in axes.ravel()[len(zs):]:
        ax.set_axis_off()

    fig.suptitle(title, color="w", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, facecolor="k")
    plt.close(fig)


def coverage(vol, zooms, title, out_path):
    """Sagittal and coronal views -- is the slab where you need it?"""
    aspect = zooms[2] / zooms[0]
    kw = display_limits(vol, phase=False)
    nx, ny = vol.shape[0], vol.shape[1]

    views = [(vol[nx // 2, :, :], "sagittal (mid)"),
             (vol[max(nx // 2 - 45, 0), :, :], "sagittal (off-mid)"),
             (vol[:, ny // 2, :], "coronal (mid)")]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5), facecolor="k")
    for ax, (sl, label) in zip(axes, views):
        ax.imshow(np.rot90(sl), aspect=aspect, **kw)
        ax.set_title(label, color="w", fontsize=10)
        ax.set_axis_off()

    fig.suptitle(title, color="w", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, facecolor="k")
    plt.close(fig)


def render_subject(sub, ses, out_root):
    """Write the four QC images. Returns a mid-brain phase slice for the
    contact sheet, or None if the subject has no usable MEGRE data."""
    session = session_dir(sub, ses)
    mag, phase = echo_files(session, "mag"), echo_files(session, "phase")

    if not mag or not phase:
        print("  sub-%s: no MEGRE mag/phase pair -- skipped" % sub)
        return None

    out_dir = os.path.join(out_root, "sub-" + sub, ses)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.join(out_dir, "sub-%s_%s_" % (sub, ses))

    first, last = min(mag), max(mag)
    last_phase = max(phase)

    def te(path):
        value = read_te_ms(path)
        return "TE %.2f ms" % value if value is not None else "TE unknown"

    m_first = nib.load(mag[first])
    m_first_data = m_first.get_fdata()

    montage(m_first_data,
            "sub-%s %s  magnitude  echo %d  (%s)" % (sub, ses, first, te(mag[first])),
            stem + "mag-echo%d.png" % first)

    m_last = nib.load(mag[last]).get_fdata()
    montage(m_last,
            "sub-%s %s  magnitude  echo %d  (%s)" % (sub, ses, last, te(mag[last])),
            stem + "mag-echo%d.png" % last)

    p_last = nib.load(phase[last_phase]).get_fdata()
    montage(p_last,
            "sub-%s %s  PHASE  echo %d  (%s)  -- look for open-ended fringe lines"
            % (sub, ses, last_phase, te(phase[last_phase])),
            stem + "phase-echo%d.png" % last_phase, phase=True)

    coverage(m_first_data, m_first.header.get_zooms()[:3],
             "sub-%s %s  magnitude echo %d -- slab coverage" % (sub, ses, first),
             stem + "coverage.png")

    print("  sub-%s: 4 images -> %s" % (sub, out_dir))

    lo, hi = brain_z_range(m_first_data)
    return p_last[:, :, (lo + hi) // 2]


def contact_sheet(slices, ses, out_path):
    """One mid-brain phase slice per subject, side by side.

    Cross-subject review is where outliers show up: the one subject whose
    phase is grainy, or wrapped far more densely than its peers.
    """
    if not slices:
        return
    cols = min(4, len(slices))
    rows = (len(slices) + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(3.4 * cols, 4.2 * rows),
                             facecolor="k", squeeze=False)
    for ax, (sub, sl) in zip(axes.ravel(), sorted(slices.items())):
        lim = float(np.abs(sl).max()) or 1.0
        ax.imshow(np.rot90(sl), cmap="gray", vmin=-lim, vmax=lim)
        ax.set_title("sub-%s" % sub, color="w", fontsize=10)
        ax.set_axis_off()
    for ax in axes.ravel()[len(slices):]:
        ax.set_axis_off()

    fig.suptitle("MEGRE phase, last echo, mid-brain -- %s" % ses,
                 color="w", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI, facecolor="k")
    plt.close(fig)
    print("\ncontact sheet (%d subjects) -> %s" % (len(slices), out_path))


def write_derivative_description(out_root):
    """BIDS requires derivatives to declare themselves."""
    path = os.path.join(out_root, "dataset_description.json")
    if os.path.exists(path):
        return
    with open(path, "w") as fh:
        json.dump({
            "Name": "BD2 MEGRE/QSM visual QC",
            "BIDSVersion": "1.8.0",
            "DatasetType": "derivative",
            "GeneratedBy": [{"Name": "7_qc_images.py"}],
        }, fh, indent=4)
        fh.write("\n")


def main():
    argv = sys.argv[1:]
    out_root = QC_DIR
    if "--outdir" in argv:
        i = argv.index("--outdir")
        if i + 1 >= len(argv):
            sys.exit("--outdir needs a path")
        out_root = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]

    if not argv:
        sys.exit("usage: %s SES [SUB ...] [--outdir DIR]"
                 % os.path.basename(sys.argv[0]))

    ses, subs = argv[0], argv[1:]
    if not subs:
        subs = discover_subjects(ses)
        if not subs:
            sys.exit("no subjects with session %r under %s" % (ses, BIDS_DIR))

    os.makedirs(out_root, exist_ok=True)
    write_derivative_description(out_root)

    print("QC images for %s (%d subject(s))" % (ses, len(subs)))
    slices = {}
    for sub in subs:
        try:
            mid = render_subject(sub, ses, out_root)
        except Exception as exc:
            print("  sub-%s: FAILED -- %s: %s" % (sub, type(exc).__name__, exc))
            continue
        if mid is not None:
            slices[sub] = mid

    contact_sheet(slices, ses, os.path.join(out_root, "phase_contact_sheet_%s.png" % ses))


main()
