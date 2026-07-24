"""Genome-browser coverage tracks — a Python port of ArchR ``plotBrowserTrack``.

Ports the coverage math from GreenleafLab/ArchR (``R/ArchRBrowser.R``:
``.regionSumArrows`` / ``.groupRegionSumArrows`` / ``.bulkTracks``):

1. Tile the region into ``tile_size`` bp bins.
2. For every cell, mark the tiles that contain a fragment **start** — binary
   per (cell, tile), exactly as ArchR does (it builds the sparse matrix from
   fragment starts and caps values at 1).
3. Sum the binary matrix over the cells of each group -> ``mat[tile, group]``
   = number of cells in the group with an insertion in that tile.
4. Scale each group by ``1e4 / group_norm_factor`` where the norm factor is the
   summed per-cell metric over the group (``nFrags`` by default; ArchR's default
   is ``ReadsInTSS``, unavailable here, so ``nFrags`` is the faithful fallback).
5. Clip to the 0.999 quantile (ArchR's default ``ylim``).

Coverage is read from 10x ``fragments.tsv.gz`` (tabix, via ``pysam``) instead of
ArchR ``.arrow`` files. Group filtering (``min_cells``) and subsampling
(``max_cells``) match ArchR's defaults.
"""

from __future__ import annotations

import glob
import re
from collections import defaultdict
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

__all__ = ["discover_fragment_files", "group_coverage", "plot_browser_track"]


def discover_fragment_files(
    raw_dir: str | Path,
    pattern: str = "*_atac_fragments.tsv.gz",
    sample_regex: str = r"_(eb\d+)_atac_fragments",
) -> dict[str, str]:
    """Map ``sample -> fragments.tsv.gz`` path from a directory of tabix files."""
    out: dict[str, str] = {}
    for p in sorted(glob.glob(str(Path(raw_dir) / pattern))):
        m = re.search(sample_regex, p)
        if m and Path(p + ".tbi").exists():
            out[m.group(1)] = p
    return out


def _default_splitter(cell: str) -> tuple[str, str]:
    """``'gse245998:eb01#AAAC...-1' -> ('eb01', 'AAAC...-1')``."""
    sample, barcode = cell.split(":", 1)[1].split("#", 1)
    return sample, barcode


def group_coverage(
    fragment_files: Mapping[str, str],
    groups: pd.Series,
    nfrags: pd.Series,
    chrom: str,
    start: int,
    end: int,
    *,
    tile_size: int = 250,
    min_cells: int = 25,
    max_cells: int = 500,
    norm: str = "nfrags",
    splitter: Callable[[str], tuple[str, str]] = _default_splitter,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray, list, float]:
    """Per-group pseudobulk coverage over ``chrom:start-end`` (ArchR port).

    Returns ``(tile_starts, mat[n_tiles, n_groups], group_names, ylim)``.
    """
    import pysam

    rng = np.random.default_rng(seed)
    groups = groups.dropna().astype(object)
    counts = groups.value_counts()
    groups = groups[groups.isin(counts[counts >= min_cells].index)]

    # subsample to max_cells per group (ArchR)
    keep: list = []
    for _, sub in groups.groupby(groups, observed=True):
        cells = sub.index.to_numpy()
        if len(cells) > max_cells:
            cells = rng.choice(cells, max_cells, replace=False)
        keep.extend(cells.tolist())
    groups = groups.loc[keep]
    uniq = sorted(groups.unique(), key=lambda x: (len(str(x)), str(x)))  # ~mixedsort
    gidx = {g: i for i, g in enumerate(uniq)}

    tile0 = start // tile_size
    n_tiles = (end // tile_size) - tile0 + 2
    tile_starts = (np.arange(n_tiles) + tile0) * tile_size

    by_sample: dict[str, dict[str, str]] = defaultdict(dict)
    for cell in groups.index:
        samp, bc = splitter(cell)
        by_sample[samp][bc] = cell

    mat = np.zeros((n_tiles, len(uniq)), dtype=np.float64)
    for samp, bc2cell in by_sample.items():
        path = fragment_files.get(samp)
        if path is None:
            continue
        tf = pysam.TabixFile(path)
        if chrom not in tf.contigs:
            continue
        cell_tiles: dict[str, set] = defaultdict(set)  # binary per (cell, tile)
        for row in tf.fetch(chrom, max(start, 0), end, parser=pysam.asTuple()):
            cell = bc2cell.get(row[3])
            if cell is None:
                continue
            ti = int(row[1]) // tile_size - tile0  # ArchR keys on the fragment start
            if 0 <= ti < n_tiles:
                cell_tiles[cell].add(ti)
        for cell, tiles in cell_tiles.items():
            j = gidx[groups.loc[cell]]
            for ti in tiles:
                mat[ti, j] += 1.0
        tf.close()

    # normalization: 1e4 / group_norm_factor
    if norm == "nfrags":
        gnf = nfrags.reindex(groups.index).groupby(groups, observed=True).sum().reindex(uniq)
    elif norm == "ncells":
        gnf = groups.value_counts().reindex(uniq)
    else:
        raise ValueError("norm must be 'nfrags' or 'ncells'")
    mat *= (1e4 / gnf.to_numpy())[np.newaxis, :]

    ylim = float(np.quantile(mat[mat > 0], 0.999)) if (mat > 0).any() else 1.0
    mat = np.clip(mat, 0.0, ylim)
    return tile_starts, mat, uniq, ylim


def plot_browser_track(
    fragment_files: Mapping[str, str],
    groups: pd.Series,
    nfrags: pd.Series,
    region: tuple[str, int, int],
    *,
    peaks: pd.DataFrame | None = None,
    genes: Mapping[str, tuple[str, int, int]] | None = None,
    highlight: tuple[int, int] | None = None,
    tile_size: int = 250,
    max_cells: int = 500,
    norm: str = "nfrags",
    palette: str = "tab20",
    group_colors: Mapping | None = None,
    order: Sequence | None = None,
    title: str = "",
    fontsize: int = 7,
    fig=None,
):
    """Draw ArchR-style tracks: per-group coverage + peaks + gene bodies.

    ``groups``/``nfrags`` are per-cell Series restricted to cells that have
    fragments (i.e. the paired GSE cells). ``peaks`` is an ATAC ``var`` frame
    with ``chrom/start/end``. ``genes`` maps name -> (chrom, start, end).
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    chrom, start, end = region
    ts, mat, uniq, ylim = group_coverage(
        fragment_files, groups, nfrags, chrom, start, end,
        tile_size=tile_size, max_cells=max_cells, norm=norm,
    )
    if order is not None:
        keep = [g for g in order if g in uniq]
        cols_idx = [uniq.index(g) for g in keep]
        mat, uniq = mat[:, cols_idx], keep

    ng = len(uniq)
    if group_colors is not None:  # match another plot's category colors (e.g. the UMAP)
        colors = [group_colors.get(g, "0.6") for g in uniq]
    else:
        colors = plt.get_cmap(palette)(np.linspace(0, 1, max(ng, 2)))
    gridspec_kw = {"height_ratios": [1] * ng + [0.3, 0.55], "hspace": 0.0}
    if fig is None:
        fig, ax = plt.subplots(
            ng + 2, 1, figsize=(9, 0.32 * ng + 2.2), sharex=True, gridspec_kw=gridspec_kw,
        )
    else:  # draw into a provided Figure / SubFigure (for composite panels)
        ax = fig.subplots(ng + 2, 1, sharex=True, gridspec_kw=gridspec_kw)
    for a, g, c in zip(ax[:ng], uniq, colors):
        if highlight:
            a.axvspan(*highlight, color="0.9", zorder=0)
        a.fill_between(ts, mat[:, uniq.index(g)], step="mid", color=c, lw=0)
        a.set_ylim(0, ylim)
        a.set_yticks([])
        a.set_ylabel(str(g), rotation=0, ha="right", va="center", fontsize=fontsize)
        for k in ("top", "right", "left"):
            a.spines[k].set_visible(False)
    ax[0].set_title(
        f"{chrom}:{start:,}-{end:,}  {title}",
        fontsize=fontsize, loc="left",
    )

    pa = ax[ng]
    pa.set_ylim(0, 1); pa.set_yticks([])
    pa.set_ylabel("Peaks", rotation=0, ha="right", va="center", fontsize=fontsize)
    if peaks is not None:
        cc = "chrom" if "chrom" in peaks.columns else peaks.columns[0]
        m = ((peaks[cc].astype(str) == chrom)
             & (peaks["end"].astype(int) > start) & (peaks["start"].astype(int) < end))
        for s, e in zip(peaks.loc[m, "start"].astype(int), peaks.loc[m, "end"].astype(int)):
            pa.add_patch(Rectangle((s, 0.25), max(e - s, tile_size), 0.5, color="#d62728", lw=0))
    for k in ("top", "right", "left"):
        pa.spines[k].set_visible(False)

    ga = ax[ng + 1]
    ga.set_ylim(0, 1); ga.set_yticks([])
    ga.set_ylabel("Genes", rotation=0, ha="right", va="center", fontsize=fontsize)
    for name, (gc, gs, ge) in (genes or {}).items():
        if gc == chrom and ge > start and gs < end:
            ga.add_patch(Rectangle((gs, 0.4), max(ge - gs, tile_size), 0.2, color="k"))
            ga.annotate(name, ((gs + ge) / 2, 0.975), ha="center", va="top",
                        fontsize=fontsize/1.5, fontstyle="italic")
    for k in ("top", "right", "left"):
        ga.spines[k].set_visible(False)
    ga.set_xlim(start, end)
    ga.set_xlabel(chrom)
    return fig, ax
