#!/usr/bin/env python3
"""Give every met entry of an adjoint ExtData.rc the time offset of the forward ExtData.rc.

Since the MAPL step_reverse fix (MAPL commit 924484fc, branch adjoint_fixes) ExtData is
evaluated at the model's own time in reverse runs, so the adjoint must read met with the
FORWARD offsets. Before that fix every adjoint met offset was the forward one minus one
dynamic step (e.g. F0;003000 -> F0;000000, PS1 0 -> 0;-003000); with the fixed binary
those shifted offsets put the adjoint met one step off.

Only entries whose file path contains MetDir are changed; everything else (restart lines,
emissions, scale factors) is left as it is.

usage: set_adjoint_met_offsets.py <forward ExtData.rc> <adjoint ExtData.rc> [--check]
       --check: change nothing, exit 1 if any met offset differs from the forward one
"""
import re
import sys


def met_offsets(path):
    off = {}
    for line in open(path):
        t = line.split('#')[0].split()
        if len(t) >= 8 and 'MetDir' in t[-1]:
            off[t[0]] = t[4]
    return off


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if len(args) != 2:
        sys.exit(__doc__)
    fwd, adj = args
    check = '--check' in sys.argv
    off = met_offsets(fwd)
    out, changed = [], []
    for line in open(adj):
        t = line.split('#')[0].split()
        if len(t) >= 8 and 'MetDir' in t[-1] and t[0] in off and t[4] != off[t[0]]:
            parts = re.split(r'(\s+)', line)
            idx = [i for i, x in enumerate(parts) if x.strip()]
            assert parts[idx[4]] == t[4], line
            changed.append(f'{t[0]}: {t[4]} -> {off[t[0]]}')
            parts[idx[4]] = off[t[0]]
            line = ''.join(parts)
        out.append(line)
    missing = sorted(set(off) - set(met_offsets(adj)))
    if missing:
        print(f'WARNING: met entries of {fwd} not in {adj}: {missing}')
    if check:
        print(f'{adj}: {len(changed)} met offsets differ from the forward ones')
        sys.exit(1 if changed else 0)
    open(adj, 'w').writelines(out)
    print(f'{adj}: {len(changed)} met offsets set to the forward value, {len(off)} met entries checked')


if __name__ == '__main__':
    main()
