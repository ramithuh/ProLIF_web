"""Extract minimal, locally sourced strict regressions from a completed audit."""
import argparse
import json
from pathlib import Path


def extract(root, output):
    cases = []
    for pdb, field, wanted in [('7gmz', 'failures', 'ImplicitHBAcceptor'),
                               ('2wyi', 'failures', 'ImplicitHBDonor'),
                               ('1wnu', 'errors', 'ImplicitHBAcceptor')]:
        native = json.loads((root / pdb / 'native.json').read_text())
        comparison = json.loads((root / pdb / 'comparison.json').read_text())
        failure = next(row for row in comparison[field] if row.get('rule') == wanted and row['lane'] == 'implicit')
        pair = next(p for p in native['pairs'] if p['ligand'] == failure['ligand'] and p['protein'] == failure['protein'])
        cases.append({'pdb': pdb, 'rule': wanted, 'kind': field,
                      'ligand': native['components'][pair['ligand']],
                      'protein': native['components'][pair['protein']],
                      'expected': pair['rules'][wanted]['expected'], 'versions': native['versions']})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({'source': 'seeded local Plinder audit 2026-09-04', 'cases': cases}, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    extract(a.root, a.output)
