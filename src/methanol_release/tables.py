"""Verify all printed Zeeman cells against the complete L26 catalogue."""
import csv,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def label(s):
    m=re.fullmatch(r'(\d+)([a-z]?)',s)
    return '$'+m[1]+('_{'+m[2]+'}' if m[2] else '')+'$'
def validate_zeeman_tables():
    rows=list(csv.DictReader((ROOT/'supplement/zeeman/top8_maser_components.csv').open()));checked=0
    blocks=list(re.finditer(r'\\begin\{table\}.*?\\end\{table\}',(ROOT/'tables/zeeman_tables.tex').read_text(),re.S))
    if len(blocks)!=27:raise ValueError('Expected 27 supplementary tables')
    groups=[]
    for b in blocks:
        family=re.search(r'\\label\{tab:sm-zeeman-([^}]+)\}',b.group())[1].replace('p','.')
        groups.append((family,b.group()))
    groups.append(('6.7',(ROOT/'tables/table-6.tex').read_text()))
    for family,block in groups:
        wanted=[r for r in rows if r['family']==family];actual=[]
        for line in block.splitlines():
            cols=[x.strip() for x in line.split('&')]
            if len(cols)==7 and re.fullmatch(r'\$\d+(?:_\{[a-z]\})?\$',cols[0]):
                cols[-1]=cols[-1].removesuffix('\\\\').strip();actual.append(cols)
        if len(actual)!=8:raise ValueError('Eight rows required')
        for cells,r in zip(actual,wanted):
            expected=[label(r['upper_label']),label(r['lower_label']),f"{float(r['delta_nu_kHz']):+.3f}",
                f"{float(r['z_Hz_per_mG']):+.4f}", f"{float(r['g1_upper']):+.5f}",
                f"{float(r['g2_lower']):+.5f}",f"{float(r['log10_A_s_inverse']):.5f}"]
            if cells!=expected:raise ValueError('Printed cell mismatch for '+family)
            checked+=1
    return dict(passed=True,checked_rows=checked,checked_cells=checked*7)
