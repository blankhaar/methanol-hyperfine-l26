"""Recalculate paper table cells, retaining their publication layout and precision."""
from pathlib import Path
import csv
import json
import re

import numpy as np

from methanol_release.forward import DATA, reduced_transform
from methanol_release.operators import THEORY, SR_TERMS, frame_parameters
from methanol_hyperfine.tensors import RHO
from methanol_line_shift.profile_estimators import broad_limit, fitted_center
from methanol_line_shift.benchmarks import BENCHMARKS

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = ('aa', 'bb', 'cc', 'ab', 'ac', 'bc')


def read_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def coefficients(problem, coordinates, theory=False):
    """Ordered M, total d and DD scales map affinely to the 98 printed quantities."""
    parameters = problem.d.default.copy() if theory else problem.unpack(coordinates)
    fractions = parameters * problem.d.coordinate_to_fraction
    m = {g: {c: np.asarray(v) * (1 + fractions[5 + SR_TERMS.index((g, c))])
             for c, v in group.items()} for g, group in THEORY['M'].items()}
    total = {g: {a: np.asarray(v, float).copy() for a, v in group.items()}
             for g, group in THEORY['total_d'].items()}
    for i, (g, axis, n) in enumerate([('CH3','a',0), ('CH3','b',0),
                                    ('OH','a',0), ('OH','b',0), ('CH3','a',1)]):
        total[g][axis][n] = parameters[i]
    lam, factor = frame_parameters()
    dp = {g: {a: v - sum(lam[j] * m[g][a+b] for j, b in enumerate('abc'))
              for a, v in group.items()} for g, group in total.items()}
    values = {}
    for group, pair, scale in [('CH3', '(1, 2)', fractions[-2]),
                               ('OH', '(0, 3)', fractions[-1])]:
        for component in COMPONENTS:
            i, j = map('abc'.index, component)
            key = 'sine_coefficients_khz' if component in ('ac','bc') else 'cosine_coefficients_khz'
            for n in range(3):
                values[f'DD_{group}_{component}_{n}'] = (1+scale) * THEORY['DD'][pair][key][n][i][j]
    for group in ('CH3', 'OH'):
        for component in COMPONENTS:
            for n in range(3):
                values[f'M_{group}_{component}_{n}'] = .5 * (m[group][component][n] + m[group][component[::-1]][n])
    for group in ('CH3', 'OH'):
        for axis in 'abc':
            for n in range(3):
                values[f'dprime_{group}_{axis}_{n}'] = dp[group][axis][n]
    values.update(A_plus=-(m['CH3']['bb'][2]+m['CH3']['cc'][2])/2,
                  A_minus=(m['CH3']['cc'][2]-m['CH3']['bb'][2])/2,
                  A_cross=(m['CH3']['bc'][2]+m['CH3']['cb'][2])/2,
                  OH_Cza=factor*dp['OH']['a'][0], CH3_Cza=factor*dp['CH3']['a'][0],
                  OH_Czz=-m['OH']['aa'][0]+(m['OH']['bb'][0]+m['OH']['cc'][0])/2+RHO*factor*dp['OH']['a'][0],
                  OH_Cxx_minus_yy=(m['OH']['cc'][0]-m['OH']['bb'][0])/2,
                  OH_Cxx_plus_yy_plus_zz=-(m['OH']['bb'][0]+m['OH']['cc'][0])/2)
    return {k: float(v) for k, v in values.items()}


def coefficient_intervals(problem, coordinates, draw_file):
    values = coefficients(problem, coordinates)
    labels = list(values)
    center = np.array(list(values.values()))
    transform = reduced_transform()[:, :13]
    derivative = np.column_stack([np.array(list(coefficients(problem, coordinates+column).values()))-center
                                  for column in transform.T])
    with np.load(draw_file, allow_pickle=False) as draws:
        if list(draws['labels']) != labels or not draws['valid'].all():
            raise ValueError('Draw labels or validity differ from the complete L26 sample')
        if np.max(abs(draws['centers']-center)) > 1e-10:
            raise ValueError('Conditional draws belong to a different molecular endpoint')
        if np.max(abs(draws['derivatives']-derivative)) > 1e-8:
            raise ValueError('Conditional coefficient derivatives disagree with the tensor map')
        samples = center + draws['coordinate_displacements'][:, :13] @ derivative.T
        nominal = samples[draws['table_sigma'] == 1]
        if len(nominal) != 512 or len(samples) != 768:
            raise ValueError('Expected 512 nominal and 256 noise-control draws')
        q16, q50, q84 = np.percentile(nominal, [16, 50, 84], axis=0)
    return values, coefficients(problem, coordinates, theory=True), {
        key: dict(center=values[key], p16=float(a), p50=float(b), p84=float(c))
        for key, a, b, c in zip(labels, q16, q50, q84)}


def rounded_like(value, example):
    """Keep the paper's chosen display precision, never its calculated value."""
    if abs(value) < 1e-10:
        value = 0.0
    if 'e' in example.lower():
        return f'{value:.6g}'
    decimals = len(example.split('.')[1]) if '.' in example else 0
    sign = '+' if example.startswith('+') else ''
    return format(value, f'{sign}.{decimals}f')


def numeric_cell(cell, value, interval=None):
    pattern = r'\$\{([^}]+)\}_\{-([^}]+)\}\^\{\+([^}]+)\}\$'
    match = re.search(pattern, cell)
    if match:
        if interval is None:
            raise ValueError('Missing uncertainty for a printed interval')
        down, up = value-interval['p16'], interval['p84']-value
        if min(down, up) < -1e-10:
            raise ValueError('Central value outside the displayed percentile interval')
        replacement = '${'+rounded_like(value, match[1])+'}_{-'+rounded_like(down, match[2])+'}^{+'+rounded_like(up, match[3])+'}$'
        return cell[:match.start()]+replacement+cell[match.end():]
    match = re.search(r'\$([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\$', cell)
    if not match:
        raise ValueError('Unrecognized numeric cell: '+cell)
    if r'\textsuperscript{b}' in cell and abs(value) > 1e-10:
        raise ValueError('A symmetry-zero table cell is no longer zero')
    return cell[:match.start(1)]+rounded_like(value, match[1])+cell[match.end(1):]


def coupling_table(text, family, values, theory, intervals):
    group, component = 'CH3', None
    output = []
    for line in text.splitlines(keepends=True):
        if r'\multicolumn{4}{l}{OH' in line or 'Methyl--OH:' in line:
            group = 'OH'
        cells = line.split('&')
        if len(cells) == 4 and re.fullmatch(r'\s*[012](?:,[012])*\s*', cells[1]):
            match = re.search(r'_\{([abc]{2})\}', cells[0])
            if match:
                component = match[1]
            if component is None:
                raise ValueError('Missing tensor component')
            keys = [f'{family}_{group}_{component}_{n.strip()}' for n in cells[1].split(',')]
            replacements = [(numeric_cell(cells[2], values[k], intervals[k]), numeric_cell(cells[3], theory[k])) for k in keys]
            if any(v != replacements[0] for v in replacements):
                raise ValueError('Grouped harmonics no longer have identical printed cells')
            cells[2], cells[3] = replacements[0]
            line = '&'.join(cells)
        output.append(line)
    return ''.join(output)


def dprime_table(text, values, intervals=None):
    group = 'CH3'
    output = []
    for line in text.splitlines(keepends=True):
        if r'\multicolumn{4}{l}{OH' in line:
            group = 'OH'
        cells = line.split('&')
        match = re.search(r"\$d'_\{([abc])\}\$", cells[0])
        if len(cells) == 4 and match:
            for n in range(3):
                key = f'dprime_{group}_{match[1]}_{n}'
                cells[n+1] = numeric_cell(cells[n+1], values[key], intervals[key] if intervals else None)
            line = '&'.join(cells)
        output.append(line)
    return ''.join(output)


def projection_table(text, keys, values, intervals):
    output, i = [], 0
    for line in text.splitlines(keepends=True):
        cells = line.split('&')
        if len(cells) == 3 and cells[0].strip().startswith(('$A_', '$C_')) and '${' in cells[2]:
            key = keys[i]; i += 1
            cells[2] = numeric_cell(cells[2], values[key], intervals[key])
            line = '&'.join(cells)
        output.append(line)
    if i != len(keys):
        raise ValueError('Projection table row count changed')
    return ''.join(output)


def zeeman_table(text, components, main=False):
    family, selected, total_rows = ('6.7' if main else None), [], 0
    def block(match):
        nonlocal family, selected, total_rows
        content = match[0]
        if not main:
            family = re.search(r'\\label\{tab:sm-zeeman-([^}]+)\}', content)[1].replace('p', '.')
        rows = {(r['upper_label'],r['lower_label']): r for r in components if r['family']==family}
        selected = []
        def row(line_match):
            nonlocal total_rows
            cells = line_match[0].split('&')
            ids = [re.sub(r'[^0-9a-z]', '', c) for c in cells[:2]]
            r = rows[tuple(ids)]; selected.append(r); total_rows += 1
            for i, (key, digits) in enumerate([('delta_nu_kHz',3),('z_Hz_per_mG',4),('g1_upper',5),('g2_lower',5),('log10_A_s_inverse',5)],2):
                cells[i] = re.sub(r'[+-]?\d+\.\d+', format(float(r[key]), f'+.{digits}f' if i<6 else f'.{digits}f'), cells[i], count=1)
            return '&'.join(cells)
        content = re.sub(r'^\$\d+(?:_\{[a-z]\})?\$\s*&.*$', row, content, flags=re.M)
        if len(selected)!=8:
            raise ValueError('Expected eight components in each magnetic table')
        ids = {(r['upper_label'],r['lower_label']) for r in selected}
        if min(float(r['A_s_inverse']) for r in selected) < max((float(r['A_s_inverse']) for k,r in rows.items() if k not in ids), default=0)*(1-1e-10):
            raise ValueError('Printed selection no longer contains the eight largest A coefficients')
        percent = 100*sum(float(r['statistical_line_fraction']) for r in selected)
        content = re.sub(r'(containing |contain\s+\$?|carry\s+\$?)(\d+\.\d+)(\\%)', lambda m:m[1]+rounded_like(percent,m[2])+m[3], content)
        return content
    pattern = r'\\begin\{table\*?\}.*?\\end\{table\*?\}'
    rendered = re.sub(pattern, block, text, flags=re.S)
    return rendered, total_rows


def response_arrays(components, transition):
    rows = [r for r in components if r['transition']==transition]
    if not rows:
        raise ValueError('Missing transition '+transition)
    f, s = np.array([[float(r['offset_kHz']),float(r['strength'])] for r in rows]).T
    return f, s/s.sum()


def line_shift_table(text, components):
    # State identities and laboratory nominal frequencies are retained inputs.
    transitions = BENCHMARKS
    output = []
    for line in text.splitlines(keepends=True):
        cells = line.split('&')
        if len(cells)==5 and re.fullmatch(r'\d+\.\d+\s*',cells[0]):
            freq = float(cells[0])*1000
            targets = [r for r in transitions if abs(r['frequency_MHz']-freq)<.001]
            if len(targets)!=1:
                raise ValueError('Ambiguous mass-ratio transition frequency')
            t=targets[0];f,s=response_arrays(components,t['id'])
            width=t['frequency_MHz']*1000/299792.458
            cells[3]=numeric_cell(cells[3],1000*broad_limit(f,s,width,'Gaussian'))
            cells[4]=numeric_cell(cells[4],1000*(fitted_center(f,s,width,'Gaussian')['center']-s@f))
            line='&'.join(cells)
        output.append(line)
    return ''.join(output)


def family_table(text, families):
    output, i = [], 0
    for line in text.splitlines(keepends=True):
        cells=line.split('&')
        if len(cells)==5 and r'\hspace{0.6em}' in cells[0]:
            r=families[i];i+=1
            count=r['N_min'] if r['N_min']==r['N_max'] else r['N_min']+'--'+r['N_max']
            value=f"{count}; {float(r['gap_min_kHz']):.1f}--{float(r['gap_max_kHz']):.1f}"
            cells[4]=re.sub(r'\d+(?:--\d+)?; [\d.]+--[\d.]+',value,cells[4],count=1)
            line='&'.join(cells)
        output.append(line)
    if i!=14:
        raise ValueError('Expected fourteen experiment families')
    return ''.join(output)


def build_tables(output, problem, coordinates, draws, numerical):
    output.mkdir()
    values,theory,intervals=coefficient_intervals(problem,coordinates,draws)
    sources={p.name:p.read_text() for p in (ROOT/'tables').glob('*.tex')}
    result={}
    for name,family in [('DD.tex','DD'),('M.tex','M')]:
        result[name]=coupling_table(sources[name],family,values,theory,intervals)
    result['dprime.tex']=dprime_table(sources['dprime.tex'],values,intervals)
    result['dprime-theory.tex']=dprime_table(sources['dprime-theory.tex'],theory)
    result['belov.tex']=projection_table(sources['belov.tex'],['A_plus','A_minus','A_cross'],values,intervals)
    result['xu.tex']=projection_table(sources['xu.tex'],['OH_Cza','CH3_Cza','OH_Czz','OH_Cxx_minus_yy','OH_Cxx_plus_yy_plus_zz'],values,intervals)
    magnetic=read_csv(numerical/'zeeman/components.csv')
    result['table-6.tex'],a=zeeman_table(sources['table-6.tex'],magnetic,True)
    result['zeeman_tables.tex'],b=zeeman_table(sources['zeeman_tables.tex'],magnetic)
    result['table-7.tex']=line_shift_table(sources['table-7.tex'],read_csv(numerical/'responses/components.csv'))
    result['table-8.tex']=family_table(sources['table-8.tex'],read_csv(numerical/'families/families.csv'))
    checks={}
    for name,text in result.items():
        (output/name).write_text(text)
        checks[name]=dict(identical_to_publication=text==sources[name])
    write_csv(output/'coupling-values-and-intervals.csv',[
        dict(parameter=k,fit_khz=v,theory_khz=theory[k],**intervals[k]) for k,v in values.items()])
    report=dict(passed=all(v['identical_to_publication'] for v in checks.values()),
                main_tables=8,supplementary_tables=27,additional_theory_table=1,
                magnetic_rows=a+b,coefficient_quantities=len(values),tables=checks,
                scope='All model/theory cells recalculated; publication captions, row selections and display precision retained; literature input cells unchanged')
    (output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
