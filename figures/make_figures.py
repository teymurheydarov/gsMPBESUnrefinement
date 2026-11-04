r"""The figures of the paper from the inputs and the runs of this repository (nothing is recomputed).

Sources:
  inputs   data/geometries/<name>.xml
  results  runs/<run>/out_final.xml; the final number of MPBES functions in the file name
           (<stem>_dof<N>) is checked against the last line of runs/<run>/out.log
Both are exported with export_visualization_mesh (examples/) into output/mesh/<name>/mesh*.txt and
drawn by
  viz      viz_main.py (the drawing script of the research code), with axes
  plain    paper_panel.py (plain line drawing; the triangle has several elements on level 0, on
           which viz_main.py does not terminate)
Figure 8 (rotor_subdomains_*) is drawn by plot_subdomains.py from the files in data/.

Output: figures/output/<stem>.png, <stem> = the file name in the paper.

usage: py figures\make_figures.py --exe <build>\bin\Release\export_visualization_mesh.exe
                                  [--path <folder with DLLs>] [stem ...]
"""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNS = ROOT / 'runs'
INPUTS = ROOT / 'data' / 'geometries'
OUT = HERE / 'output'

# stem -> (source, drawer); source: 'run:<run in runs/>' or 'input:<xml name in data/geometries>'
FIGURES = {
    # inputs
    'mask_bands_original':                         ('input:mask_bands_L3', 'viz'),
    'yeti_mp2_thb_approximation_fine_L3_original': ('input:yeti_mp2_thb_approximation_fine_L3_topology', 'viz'),
    'joystick_approximation_fine_L3_NLO_original': ('input:joystick_approximation_fine_L3_NLO', 'viz'),
    'rotorRE_original':                            ('input:rotorRE', 'viz'),
    'tri_wedge_original':                          ('input:triangle_wedge_L5', 'plain'),
    # results
    'yeti_eg001_ef001_dof239':      ('run:yeti_eg001_ef001', 'viz'),
    'yeti_eg003_ef003_dof157':      ('run:yeti_eg003_ef003', 'viz'),
    'yeti_eg006_ef006_dof122':      ('run:yeti_eg006_ef006', 'viz'),
    'yeti_F_ef010_dof178':          ('run:yeti_F_ef010', 'viz'),
    'yeti_F_ef030_dof128':          ('run:yeti_F_ef030', 'viz'),
    'yeti_F_ef060_dof117':          ('run:yeti_F_ef060', 'viz'),
    'rotor_eg0001_ef0001_dof476':   ('run:rotor_eg0001_ef0001', 'viz'),
    'rotor_nog_ef0001_dof242':      ('run:rotor_nog_ef0001', 'viz'),
    'mask_bands_ls_dof568':         ('run:mask_bands_ls_only', 'viz'),
    'mask_bands_lo_dof562':         ('run:mask_bands_lo_only', 'viz'),
    'mask_bands_nlo_full_dof351':   ('run:mask_bands_nlo_full', 'viz'),
    'tri_wedge_ls_dof1994':         ('run:triangle_wedge_ls_only', 'plain'),
    'tri_wedge_lo_dof349':          ('run:triangle_wedge_lo_only', 'plain'),
    'joystick_grendas_dof512':      ('run:joystick_grendas', 'viz'),
    'rotor_lambda0_ef00002_dof356': ('run:rotor_lambda0_ef00002', 'viz'),
    'rotor_lambda2_ef00002_dof405': ('run:rotor_lambda2_ef00002', 'viz'),
    'rotor_global_ef00002_dof478':  ('run:rotor_global_ef00002', 'viz'),
}
SUBDOMAINS = ['rotor_subdomains_' + n for n in ('geometry', 'lambda0', 'lambda2', 'global')]


def final_count(run):
    last = (RUNS / run / 'out.log').read_text(encoding='utf-8').strip().splitlines()[-1]
    m = re.search(r'-> (\d+) MPBES functions', last)
    if not m:
        sys.exit(f'ERROR: {run}: no final count in the last line of out.log: {last}')
    return int(m.group(1))


def mesh_file(source, exe, env):
    kind, name = source.split(':', 1)
    xml = RUNS / name / 'out_final.xml' if kind == 'run' else INPUTS / f'{name}.xml'
    prefix = OUT / 'mesh' / name / 'mesh'   # short: Windows limits paths to 260 characters
    if not prefix.with_suffix('.txt').exists():
        prefix.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([exe, '-f', str(xml), '-o', str(prefix)], env=env, check=True,
                       stdout=subprocess.DEVNULL)
        # the export does not report files it could not open
        longest = Path(str(prefix) + '_interface_meta.txt')
        if not longest.exists():
            sys.exit(f'ERROR: {name}: export wrote no {longest} '
                     f'({len(str(longest))} characters; Windows allows 260)')
    return prefix.with_suffix('.txt')


def draw_viz(stem, mesh, png):
    r = subprocess.run([sys.executable, str(HERE / 'viz_main.py'), str(mesh), str(png)],
                       capture_output=True, text=True)
    points = png.with_name(png.stem + '_points.png')   # point-cloud by-product of viz_main.py
    if points.exists():
        points.unlink()
    if r.returncode != 0 or not png.exists() or '[classify] all' not in r.stdout:
        sys.exit(f'ERROR: {stem}: viz_main.py failed or lines unmatched\n'
                 f'{r.stdout[-2000:]}\n{r.stderr[-2000:]}')


def draw_plain(stem, mesh, png):
    subprocess.run([sys.executable, str(HERE / 'paper_panel.py'), str(png), str(mesh)], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--exe', required=True, help='export_visualization_mesh executable')
    ap.add_argument('--path', action='append', default=[], help='folder prepended to PATH (DLLs)')
    ap.add_argument('stems', nargs='*', help='figures to draw (default: all)')
    a = ap.parse_args()
    exe = str(Path(a.exe).resolve())
    env = dict(os.environ, PATH=os.pathsep.join([str(Path(p).resolve()) for p in a.path]
                                                + [os.environ['PATH']]))
    OUT.mkdir(exist_ok=True)
    for stem, (source, drawer) in FIGURES.items():
        if a.stems and stem not in a.stems:
            continue
        if source.startswith('run:'):
            n = final_count(source[4:])
            if not stem.endswith(f'_dof{n}'):
                sys.exit(f'ERROR: {stem}: run {source[4:]} ends with {n} MPBES functions')
        png = OUT / f'{stem}.png'
        (draw_viz if drawer == 'viz' else draw_plain)(stem, mesh_file(source, exe, env), png)
        print(f'{stem}: {png}', flush=True)
    if not a.stems or any(s in a.stems for s in SUBDOMAINS):
        subprocess.run([sys.executable, str(HERE / 'plot_subdomains.py')], check=True)


if __name__ == '__main__':
    main()
