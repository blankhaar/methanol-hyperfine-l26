"""Plot Figures 2--4 at the selected L26 endpoint. Molecular parameters fixed.

The frozen numerical inputs are supplied in evidence/figure_data/.
Run from any directory with NumPy and Matplotlib installed.
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FixedLocator
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'evidence/figure_data'
OUT = ROOT / 'manuscript/figures'
OUT.mkdir(parents=True, exist_ok=True)
A = np.load(DATA / 'figure_data.npz', allow_pickle=False)
META = json.loads((DATA / 'figure_metadata.json').read_text())
PROFILES = {r['id']: r for r in META['profiles']}
INPUT_HASHES = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in DATA.iterdir() if p.is_file()}
TEAL, GRAY = '#007f89', '.48'
STACKS = []
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})


def style(ax):
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.6)
    ax.tick_params(length=3, width=.6, pad=2, direction='in', right=True)
    ax.margins(x=.015)


def stack_style(ax, rx):
    for panel in (ax, rx):
        style(panel)
    ax.tick_params(axis='x', labelbottom=False, bottom=False, top=True)
    rx.tick_params(axis='x', top=False, bottom=True)
    ax.yaxis.set_major_locator(MaxNLocator(3, prune='lower'))
    rx.yaxis.set_major_locator(MaxNLocator(3))
    STACKS.append((ax, rx))


def pair(fig, cell, name, title, legend=False):
    grid = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=0)
    ax = fig.add_subplot(grid[0])
    rx = fig.add_subplot(grid[1], sharex=ax)
    x, y, pred, mask = [A[name + '__' + k] for k in ['x', 'y', 'pred', 'mask']]
    ax.plot(x, y, color=GRAY, lw=.65, label='Experiment')
    ax.plot(x, pred, color=TEAL, lw=1.1, label='L26')
    cuts = np.flatnonzero(np.diff(np.r_[False, ~mask, False]))
    for start, stop in zip(cuts[::2], cuts[1::2]):
        for panel in (ax, rx):
            panel.axvspan(x[start], x[stop - 1], color='.87', zorder=-2)
    rx.plot(x, np.where(mask, y - pred, np.nan), color=GRAY, lw=.55)
    rx.axhline(0, color='.7', lw=.6)
    ax.set_title(title, loc='left', fontsize=8.4, pad=5)
    ax.set_ylabel('Signal')
    rx.set_ylabel('O-C')
    rx.set_xlabel('Frequency offset (kHz)', labelpad=2)
    if legend:
        ax.legend(frameon=False, fontsize=7.6, loc='upper left', handlelength=1.5)
    stack_style(ax, rx)
    return ax, rx


def save(fig, name):
    fig.canvas.draw()
    for ax, rx in STACKS:
        if ax.figure is not fig:
            continue
        # Keep y tick labels away from the shared seam; MaxNLocator pruning
        # alone can prune an off-screen tick and leave a near-boundary label.
        lo, hi = ax.get_ylim()
        ax.yaxis.set_major_locator(FixedLocator([v for v in ax.get_yticks()
            if .10 < (v-lo)/(hi-lo) < .98]))
        lo, hi = rx.get_ylim()
        rx.yaxis.set_major_locator(FixedLocator([v for v in rx.get_yticks()
            if .02 < (v-lo)/(hi-lo) < .83]))
        p, q = ax.get_position(), rx.get_position()
        assert abs(p.y0 - q.y1) < 1e-10
        assert abs(p.x0 - q.x0) < 1e-10 and abs(p.x1 - q.x1) < 1e-10
    fig.savefig(OUT / (name + '.pdf'), metadata={
        'Title': name, 'Subject': 'L26 predictions'})
    plt.close(fig)


# Figure 4: the same six selected low-frequency spectra.
fig = plt.figure(figsize=(6.69, 7.3))
grid = fig.add_gridspec(3, 2, left=.09, right=.985, bottom=.07, top=.965,
                      hspace=.5, wspace=.29)
ids = ['C2015-p01-2', 'C2015-p02-2', 'C2015-p03-2',
       'C2015-p05-2', 'C2015-p07-2', 'C2015-p10-1']
for i, name in enumerate(ids):
    g = PROFILES[name]
    ax, rx = pair(fig, grid[i // 2, i % 2], name,
                  f"({chr(97+i)}) {g['lab']}: {g['frequency_mhz']/1000:.4f} GHz {g['species']}", i == 0)
    if i == 0:
        # Requested symmetric viewing window only; shared axes leave the
        # original sample coordinates, predictions and residuals untouched.
        ax.set_xlim(-25, 25)
        assert tuple(rx.get_xlim()) == (-25., 25.)
save(fig, 'fig1-low-j')

# Figure 2: Belov separations and two Xu profiles.
rows = list(csv.DictReader((DATA / 'belov.csv').open()))
# Belov et al. (2016), Table IV, Fit A(2); source transcription checked in review.
published = json.loads((DATA / 'belov-published-table-iv.json').read_text())['rows']
by_state = {(r['ju'],r['ku'],r['jl'],r['kl']):r for r in published}
assert len(by_state) == len(rows) == 57
for r in rows:
    p = by_state[tuple(int(r[k]) for k in ['ju','ku','jl','kl'])]
    assert abs(float(r['observed_khz'])-abs(p['signed_observed_khz'])) < 1e-10
    assert (r['included'] == 'True') == bool(p['published_weight'])
    r['published_A2_khz'] = abs(p['fit_A2_khz'])
assert sum(r['included'] == 'True' for r in rows) == 55
fig = plt.figure(figsize=(6.69, 5.5))
outer = fig.add_gridspec(2, 1, left=.095, right=.985, bottom=.09, top=.955,
                       height_ratios=[1.35, 1], hspace=.48)
top = outer[0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0)
ax = fig.add_subplot(top[0])
rx = fig.add_subplot(top[1], sharex=ax)
families = [((1,-1,-2), r'R: $-1\leftarrow-2$', '#36894b', '^'),
            ((0,0,-1), r'Q: $0\leftarrow-1$', '#3b63b8', 'o'),
            ((0,2,-1), r'Q: $2\leftarrow-1$', '#d54e4e', 'o'),
            ((0,3,-2), r'Q: $3\leftarrow-2$', '#8754ac', '^'),
            ((0,-2,1), r'Q: $-2\leftarrow1$', '#40494f', 's')]
for key, label, col, marker in families:
    rr = sorted([r for r in rows if (int(r['ju'])-int(r['jl']), int(r['ku']), int(r['kl'])) == key],
                key=lambda r: int(r['jl']))
    j = np.array([int(r['jl']) for r in rr])
    y = np.array([float(r['observed_khz']) for r in rr])
    v = np.array([float(r['predicted_khz']) for r in rr])
    keep = np.array([r['included'] == 'True' for r in rr])
    ax.plot(j, v, color=col, lw=1)
    ax.plot(j, [r['published_A2_khz'] for r in rr], color=col, lw=1.15, ls=(0,(3,2)))
    ax.scatter(j[keep], y[keep], s=17, marker=marker, facecolor='white', edgecolor=col, lw=.8, label=label, zorder=3)
    ax.scatter(j[~keep], y[~keep], s=24, marker='x', color=col, lw=1)
    rx.scatter(j[keep], (y-v)[keep], s=12, marker=marker, color=col)
    rx.scatter(j[~keep], (y-v)[~keep], s=22, marker='x', color=col, lw=1)
ax.set_title('(a) Belov: high-J E doublet separations', loc='left', fontsize=8.4, pad=5)
ax.set_ylabel('Splitting (kHz)')
ax.set_ylim(23, 79)
family_legend = ax.legend(frameon=False, loc='upper left', fontsize=7.3, ncol=2, handletextpad=.4, columnspacing=.7)
ax.add_artist(family_legend)
ax.legend(handles=[Line2D([],[],color='.2',lw=1,label='L26'),
                   Line2D([],[],color='.2',lw=1.15,ls=(0,(3,2)),label='Belov Fit A(2)')],
          loc='upper right', frameon=True, facecolor='white', edgecolor='none',
          fontsize=7, handlelength=2.5, framealpha=.95)
rx.axhspan(-1, 1, color='#edf3f3', zorder=-1)  # Visual ±1-kHz guide; no uncertainty meaning.
rx.axhline(0, color='.6', lw=.6)
rx.set_ylim(-2.2, 3.3)
rx.set_ylabel('O-C (kHz)')
rx.set_xlabel('Lower-state J', labelpad=1)
rx.set_xticks(np.arange(14, 35, 2))
rx.set_xlim(12.5, 34.7)
stack_style(ax, rx)
bottom = outer[1].subgridspec(1, 2, wspace=.29)
for i, (name, title) in enumerate([
        ('Xu-F2-34', r'(b) Xu: $v_t=1$, Q(3), $K=3\leftarrow2$'),
        ('Xu-F2-48', r'(c) Xu: $v_t=1$, R(2), $K=3\leftarrow2$')]):
    pair(fig, bottom[i], name, title, legend=(i == 0))
save(fig, 'fig2-belov-xu')

# Figure 3: retain the crossover-inclusive fitted spectrum; omit the comparison curve.
fig = plt.figure(figsize=(3.37, 3.1))
grid = fig.add_gridspec(2, 1, left=.18, right=.97, bottom=.17, top=.97,
                      height_ratios=[3, 1], hspace=0)
ax = fig.add_subplot(grid[0])
rx = fig.add_subplot(grid[1], sharex=ax)
x = A['triplet__x'] + float(META['triplet_source_origin_khz'])
y, pr = A['triplet__y'], A['triplet__with']
ax.scatter(x, y, s=6, color=GRAY, label='Experiment', zorder=3)
ax.plot(x, pr, lw=1.15, color=TEAL, label='L26')
ax.set_ylabel('Signal')
ax.set_ylim(min(y.min(), pr.min())-.06, max(y.max(), pr.max())+.38)
ax.legend(loc='upper left', frameon=False, fontsize=7.4, handlelength=1.8,
          labelspacing=.3, borderaxespad=.5)
rx.scatter(x, y-pr, s=4, color=GRAY)
rx.axhline(0, color='.65', lw=.6)
rx.set_ylabel('O-C')
rx.set_xlabel('Frequency offset (kHz)', labelpad=3)
rx.set_xlim(-155, 155)
rx.set_xticks([-150, -75, 0, 75, 150])
stack_style(ax, rx)
save(fig, 'fig3-belov-triplet')

assert INPUT_HASHES == {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in DATA.iterdir() if p.is_file()}
(ROOT / 'spectral-figure-validation.json').write_text(json.dumps({
    'input_sha256': INPUT_HASHES,
    'endpoint_sha256': META['endpoint_sha256'],
    'molecular_refit': False, 'response_refit_during_plotting': False,
    'triplet_response_refitted_upstream_at_fixed_endpoint': True,
    'all_stack_gaps_zero': True, 'stack_count': len(STACKS),
    'figure2_Belov_comparison': '57 absolute Table IV Fit A(2) values, identity/observation/weight matched; 55 retained, 2 excluded',
    'figure2_shading': 'Residual interval [-1,+1] kHz, visual guide only',
    'figure3_curve': 'triplet__with', 'figure3_comparison_curve_present': False,
    'figure3_points': len(y), 'residual': 'observed minus crossover-inclusive prediction',
    'figure4a_view_xlim_khz': [-25, 25],
    'figure4a_native_sample_range_khz': [float(A['C2015-p01-2__x'].min()),
                                       float(A['C2015-p01-2__x'].max())],
    'figure4a_samples_outside_view': int(np.count_nonzero(
        (A['C2015-p01-2__x'] < -25) | (A['C2015-p01-2__x'] > 25))),
}, indent=2)+'\n')
print('Rebuilt Figures 2, 3 and 4 using L26.')
