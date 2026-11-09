"""Plot fresh L26 line-pull calculations in the original Figure 5 style."""
from pathlib import Path
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'evidence/line_pull_data'
OUT = ROOT / 'manuscript/figures'
C_KMS = 299792.458
TARGETS = [('12.2GHz', 12178.596415), ('96.741GHz', 96741.375)]
STYLES = {
    'Gaussian': ('#007F7B', '-', 'Gaussian'),
    'Lorentz': ('#C06C20', '--', 'Lorentzian'),
    'dispersion_1f': ('#96629A', '-.', '1f dispersion'),
    'absorption_2f': ('#335C9B', ':', '2f absorption'),
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs = ['line_pull_curves.csv', 'line_pull_figure_validation.json', 'components.csv']
    hashes = {name: sha256(DATA / name) for name in inputs}
    validation = json.loads((DATA / 'line_pull_figure_validation.json').read_text())
    assert hashes['components.csv'] == validation['source_catalog_sha256']
    rows = list(csv.DictReader((DATA / 'line_pull_curves.csv').open()))
    assert len(rows) == validation['curve_rows'] == 1464
    grids = {(target, response): [r for r in rows if r['target'] == target and r['response'] == response]
             for target, _ in TARGETS for response in STYLES}
    assert all(len(r) == 183 for r in grids.values())
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'axes.linewidth': .7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(1, 2, figsize=(8.9, 4.15), sharey=True)
    fig.subplots_adjust(left=.085, right=.985, bottom=.15, top=.97, wspace=.10)
    conversion_checks = []
    for i, (ax, (target, frequency_mhz)) in enumerate(zip(axes, TARGETS)):
        factor = C_KMS / (frequency_mhz * 1000)
        for response, (color, ls, label) in STYLES.items():
            selected = grids[target, response]
            widths_khz = np.array([float(r['component_absorption_FWHM_kHz']) for r in selected])
            widths_kms = widths_khz * factor
            saved_velocity = np.array([float(r['component_velocity_FWHM_kms']) for r in selected])
            assert np.allclose(widths_kms, saved_velocity, rtol=2e-15, atol=0)
            pull = np.array([float(r['fitted_center_minus_centroid_Hz']) for r in selected])
            assert np.all(pull > 0) if i == 0 else np.all(pull < 0)
            ax.plot(widths_kms, abs(pull), color=color, ls=ls, lw=1.85, label=label)
            conversion_checks.append({'target': target, 'response': response,
                'points': len(selected), 'velocity_conversion_max_abs_kms': float(np.max(abs(widths_kms-saved_velocity)))})
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(30*factor, 2000*factor)
        ax.set_ylim(.015, 5000)
        ticks = [1, 3, 10, 30] if i == 0 else [.1, .3, 1, 3]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_locator(FixedLocator([.1, 1, 10, 100, 1000]))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:g}'))
        ax.grid(axis='y', which='major', color='#DDE2E5', lw=.6)
        ax.tick_params(which='both', direction='in', top=True, right=True)
        for spine in ax.spines.values():
            spine.set_visible(True)
        ax.set_xlabel(r'$W$ (km/s)', labelpad=7)
        title = '(a) 12.2 GHz' if i == 0 else '(b) 96.741 GHz'
        ax.text(.97, .97, title, transform=ax.transAxes, ha='right', va='top', fontsize=11,
                bbox={'facecolor':'white', 'edgecolor':'none', 'pad':1.5})
    axes[0].set_ylabel(r'$|c|$ (Hz)', labelpad=8)
    axes[0].legend(loc='lower left', fontsize=8, frameon=True, facecolor='white',
                   edgecolor='none', framealpha=1, handlelength=2.6,
                   borderpad=.35, labelspacing=.35)
    teal = STYLES['Gaussian'][0]
    factor = C_KMS/(TARGETS[0][1]*1000)
    for width, xytext in [(100., (70., 3.0)), (300., (300., .3))]:
        r = next(r for r in grids['12.2GHz', 'Gaussian']
                 if float(r['component_absorption_FWHM_kHz']) == width)
        pull = float(r['fitted_center_minus_centroid_Hz'])
        axes[0].plot(width*factor, abs(pull), 'o', ms=4.4, mfc='white', mec=teal, mew=1.2, zorder=5)
        axes[0].annotate(f'{pull:+.1f} Hz', (width*factor, abs(pull)),
                         xytext=(xytext[0]*factor, xytext[1]), textcoords='data',
                         fontsize=9, color=teal,
                         arrowprops={'arrowstyle':'-', 'color':teal, 'lw':.7},
                         bbox={'facecolor':'white', 'edgecolor':'none', 'pad':1})
    cloud = validation['cloud_illustration']
    factor = C_KMS/(TARGETS[1][1]*1000)
    cloud_width = cloud['component_absorption_FWHM_kHz']*factor
    cloud_pull = cloud['fitted_center_minus_centroid_Hz']
    axes[1].plot(cloud_width, abs(cloud_pull), 'o', ms=4.8, mfc='white', mec=teal, mew=1.2, zorder=5)
    axes[1].annotate(f'{cloud_pull:+.0f} Hz', (cloud_width, abs(cloud_pull)),
                     xytext=(75*factor, 25), textcoords='data', fontsize=9, color=teal,
                     arrowprops={'arrowstyle':'-', 'color':teal, 'lw':.8}, ha='left', va='top')
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ['pdf', 'png']:
        fig.savefig(OUT/f'fig4-hyperfine-line-pulls.{ext}', dpi=240)
    plt.close(fig)
    assert hashes == {name: sha256(DATA/name) for name in inputs}
    (ROOT/'line-pull-display-validation.json').write_text(json.dumps({
        'status': 'passed', 'input_sha256': hashes, 'curve_rows': len(rows),
        'endpoint_sha256': validation['source_endpoint_sha256'],
        'molecular_refit': False, 'response_refit': False,
        'axis_x': 'individual component absorption FWHM in km/s',
        'axis_y': 'absolute fitted-center shift in Hz',
        'conversion_checks': conversion_checks,
        'cloud_example': cloud,
    }, indent=2)+'\n')
    print('Rebuilt Figure 5 from 1464 newly calculated L26 response ordinates.')


if __name__ == '__main__':
    main()
