#!/usr/bin/env python3
"""Check that the adjoint flux is the forward flux without the control variable.

dJ/dsigma = sum lambda * F_adj, where F_adj (the HEMCO flux the adjoint uses) must be
the forward flux with every factor EXCEPT sigma (HEMCO scale factor 750).  For every
base-emission line of the forward HEMCO_Config.rc that carries scale 750:
  - the adjoint must have the same container with the same file/variable,
  - its scale list must be the forward list with 750 removed,
  - every other scale factor in that list must be defined identically in both files.
A mismatch (e.g. the OSSE prior-error factor missing in the adjoint) makes the
gradient wrong without any error, so this script exits 1 on any difference.

usage: check_adjoint_flux_config.py <forward HEMCO_Config.rc> <adjoint HEMCO_Config.rc> [--sigma-id 750]
"""
import sys


def parse(path):
    base, scales = {}, {}
    for line in open(path):
        s = line.split('#')[0].strip()
        t = s.split()
        if len(t) >= 12 and t[0] == '0':           # base emission: 0 NAME file var time C xy unit spc scal cat hier
            base[t[1]] = t
        elif len(t) >= 2 and t[0].isdigit() and t[0] != '0' and not s.startswith('('):
            scales.setdefault(t[0], ' '.join(t))   # first definition wins (HEMCO order)
    return base, scales


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        sys.exit(__doc__)
    sid = sys.argv[sys.argv.index('--sigma-id') + 1] if '--sigma-id' in sys.argv else '750'
    fb, fs = parse(args[0])
    ab, as_ = parse(args[1])
    problems, n = [], 0
    for name, t in fb.items():
        flist = t[9].split('/')
        if sid not in flist:
            continue
        n += 1
        want = [x for x in flist if x != sid]
        if name not in ab:
            problems.append(f'{name}: scaled by sigma in the forward run but missing in the adjoint')
            continue
        a = ab[name]
        if a[2:4] != t[2:4]:
            problems.append(f'{name}: file/variable differ: forward {t[2:4]} adjoint {a[2:4]}')
        alist = [x for x in a[9].split('/') if x != '-']
        if alist != want:
            problems.append(f'{name}: adjoint scale list {a[9]} != forward {t[9]} without {sid} '
                            f'({"/".join(want) or "-"})')
        for x in want:
            if fs.get(x) != as_.get(x):
                problems.append(f'scale factor {x} differs:\n      forward: {fs.get(x)}\n      adjoint: {as_.get(x)}')
    if n == 0:
        problems.append(f'no forward flux carries the sigma factor {sid}')
    if problems:
        print('ADJOINT FLUX CONFIG CHECK FAILED:\n  - ' + '\n  - '.join(problems))
        sys.exit(1)
    print(f'adjoint flux config OK: {n} sigma-scaled flux(es), adjoint = forward without factor {sid}')


if __name__ == '__main__':
    main()
