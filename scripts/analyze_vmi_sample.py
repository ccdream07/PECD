"""Analyze one Andor/VMI image with PyAbel rBasex.

This is a first-pass diagnostic tool for the sample in ``test``. It does not
claim to calculate PECD from a single image: quantitative PECD requires paired
LCP/RCP images with their helicity labels. The odd Legendre contribution from
one image is reported only as a symmetry/QC diagnostic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from abel.rbasex import rbasex_transform
from abel.tools import center as abel_center


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image", type=Path, help="raw 2-D numeric image, e.g. test/200-130-60-2.dat")
    parser.add_argument("--outdir", type=Path, default=None, help="output directory")
    parser.add_argument("--center-row", type=float, default=None)
    parser.add_argument("--center-col", type=float, default=None)
    parser.add_argument("--order", type=int, default=4, help="highest Legendre order")
    parser.add_argument("--rmax", default="MIN", help="rBasex radial limit")
    parser.add_argument(
        "--vendor-aniso",
        type=Path,
        default=None,
        help="optional Basex aniso.dat for side-by-side diagnostic plotting",
    )
    return parser.parse_args()


def choose_origin(image: np.ndarray, row: float | None, col: float | None) -> tuple[float, float]:
    if row is not None and col is not None:
        return int(round(row)), int(round(col))
    if row is not None or col is not None:
        raise ValueError("--center-row and --center-col must be supplied together")
    origin = abel_center.find_origin_by_convolution(image, round_output=False)
    return int(round(float(origin[0]))), int(round(float(origin[1])))


def save_image(path: Path, image: np.ndarray, title: str, diverging: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    if diverging:
        vmax = float(np.nanpercentile(np.abs(image), 99.5))
        vmax = max(vmax, np.finfo(float).eps)
        im = ax.imshow(image, origin="upper", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    else:
        im = ax.imshow(np.log1p(np.clip(image, 0, None)), origin="upper", cmap="magma")
    ax.set_title(title)
    ax.set_xlabel("camera column")
    ax.set_ylabel("camera row")
    fig.colorbar(im, ax=ax, label="log(1 + counts)" if not diverging else "reconstructed intensity")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    image_path = args.image.resolve()
    if not image_path.exists():
        raise FileNotFoundError(image_path)
    outdir = (args.outdir or image_path.parent / f"{image_path.stem}_analysis").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    image = np.loadtxt(image_path, dtype=float)
    if image.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {image.shape}")
    image = np.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    origin = choose_origin(image, args.center_row, args.center_col)

    reconstruction, results = rbasex_transform(
        image,
        origin=origin,
        rmax=args.rmax,
        order=args.order,
        odd=True,
        out="same",
        verbose=False,
    )

    # Ibeta rows are: I(r), beta_1(r), beta_2(r), ... beta_order(r).
    ibeta = results.Ibeta(window=3)
    radius = np.asarray(results.r, dtype=float)
    intensity = np.asarray(ibeta[0], dtype=float)
    beta = np.asarray(ibeta[1:], dtype=float)

    columns = ["radius_px", "intensity"] + [f"beta_{i}" for i in range(1, beta.shape[0] + 1)]
    table = np.column_stack([radius, intensity, beta.T])
    np.savetxt(outdir / "rbasex_radial_anisotropy.csv", table, delimiter=",", header=",".join(columns), comments="")
    np.save(outdir / "rbasex_reconstruction.npy", reconstruction)

    # Very weak radii are dominated by background and inversion noise. Keep
    # the full curves, but use a conservative 1% intensity threshold for QC.
    valid_fraction = 0.01
    positive = intensity > max(float(np.nanmax(intensity)) * valid_fraction, 0.0)
    if np.any(positive):
        odd_abs = np.nanmedian(np.abs(beta[[0, 2]][:, positive])) if beta.shape[0] >= 3 else float("nan")
        even_abs = np.nanmedian(np.abs(beta[[1, 3]][:, positive])) if beta.shape[0] >= 4 else float("nan")
        odd_component = 2 * beta[0] - 0.5 * beta[2] if beta.shape[0] >= 3 else np.full_like(radius, np.nan)
    else:
        odd_abs = even_abs = float("nan")
        odd_component = np.full_like(radius, np.nan)

    np.savetxt(
        outdir / "single_image_odd_component.csv",
        np.column_stack([radius, odd_component]),
        delimiter=",",
        header="radius_px,2beta1_minus_half_beta3_not_PECD",
        comments="",
    )

    save_image(outdir / "raw_image.png", image, "Raw Andor/VMI image")
    save_image(outdir / "rbasex_reconstruction.png", reconstruction, "PyAbel rBasex reconstruction", diverging=True)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(radius, np.maximum(intensity, 0), lw=1.2)
    ax.set_xlabel("radius (pixel)")
    ax.set_ylabel("radial intensity")
    ax.set_title("Radial photoelectron distribution")
    ax.grid(alpha=0.25)
    fig.savefig(outdir / "radial_spectrum.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)
    for i in range(beta.shape[0]):
        axes[0].plot(radius, beta[i], label=rf"$\beta_{i+1}$", lw=1)
    axes[0].axhline(0, color="black", lw=0.7)
    axes[0].set_ylabel("anisotropy coefficient")
    axes[0].legend(ncol=4, fontsize=9)
    axes[0].grid(alpha=0.25)
    axes[1].plot(radius, odd_component, color="tab:red", lw=1.1)
    axes[1].axhline(0, color="black", lw=0.7)
    axes[1].set_xlabel("radius (pixel)")
    axes[1].set_ylabel(r"$2\beta_1-\beta_3/2$")
    axes[1].set_title("Single-image odd component (not PECD)")
    axes[1].grid(alpha=0.25)
    fig.savefig(outdir / "anisotropy_qc.png", dpi=160)
    plt.close(fig)

    vendor_path = args.vendor_aniso
    if vendor_path is None:
        candidate = image_path.with_name(f"{image_path.stem}_aniso.dat")
        vendor_path = candidate if candidate.exists() else None
    vendor_summary = None
    if vendor_path is not None and vendor_path.exists():
        vendor = np.loadtxt(vendor_path, dtype=float)
        vendor_summary = {"path": str(vendor_path.resolve()), "shape": list(vendor.shape), "columns": int(vendor.shape[1])}
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        for j in range(1, vendor.shape[1]):
            ax.plot(vendor[:, 0], vendor[:, j], label=f"Basex column {j}", lw=1)
        ax.set_xlabel("radius (pixel)")
        ax.set_ylabel("vendor Basex value")
        ax.set_title("Vendor Basex anisotropy output; column meanings not assumed")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.25)
        fig.savefig(outdir / "vendor_aniso.png", dpi=160)
        plt.close(fig)

    summary = {
        "input": str(image_path),
        "input_shape": list(image.shape),
        "origin_row_col": list(origin),
        "method": "PyAbel rBasex",
        "order": int(args.order),
        "rmax": args.rmax,
        "radial_points": int(radius.size),
        "qc_intensity_threshold_fraction": valid_fraction,
        "median_abs_odd_coefficients": None if not np.isfinite(odd_abs) else float(odd_abs),
        "median_abs_even_coefficients": None if not np.isfinite(even_abs) else float(even_abs),
        "interpretation": "A single image cannot provide quantitative LCP/RCP PECD. The odd component is a symmetry/QC diagnostic only.",
        "vendor_aniso": vendor_summary,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
