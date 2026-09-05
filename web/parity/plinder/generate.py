"""Native ProLIF oracle on frozen local Plinder PDB/SDF systems; no downloads."""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_case(folder):
    import numpy as np
    import prolif
    import rdkit
    from rdkit import Chem
    from prolif import interactions as pi
    from prolif.molecule import Molecule
    from prolif.residue import Residue

    names = ['Hydrophobic', 'HBAcceptor', 'HBDonor', 'ImplicitHBAcceptor',
             'ImplicitHBDonor', 'XBAcceptor', 'XBDonor', 'Cationic', 'Anionic',
             'CationPi', 'PiCation', 'PiStacking', 'MetalDonor', 'MetalAcceptor', 'VdWContact']
    rules = {name: getattr(pi, name)() for name in names}
    rules['VdWContactCSD'] = pi.VdWContact(preset='csd')
    result = {'versions': {'prolif': prolif.__version__, 'rdkit': rdkit.__version__,
                         'prolifPath': prolif.__file__, 'python': sys.version},
              'ruleSourceHashes': {name: sha(inspect.getfile(type(rule))) for name, rule in rules.items()},
              'components': {}, 'pairs': [], 'errors': [], 'lanes': ['implicit', 'explicit'],
              'rules': list(rules), 'selection': 'all ligand SDFs; receptor residues within 8 A of any ligand atom'}

    def normalize(name, meta):
        def number(x):
            if isinstance(x, dict): return {k: number(v) for k, v in x.items()}
            if isinstance(x, (tuple, list, np.ndarray)): return [number(v) for v in x]
            if isinstance(x, (np.integer,)): return int(x)
            if isinstance(x, (np.floating,)): return float(x)
            return x
        return {'interaction': 'VdWContact' if name == 'VdWContactCSD' else name,
                'indices': number(meta['indices']), 'parentIndices': number(meta['parent_indices']),
                'distance': float(meta['distance']),
                'geometry': number({k: v for k, v in meta.items() if k not in ['indices', 'parent_indices', 'distance']})}

    def component(residue, key):
        if key in result['components']: return cache[key]
        # Both engines read the exact same serialized graph, coordinates and H state.
        block = Chem.MolToMolBlock(residue)
        molecule = Chem.MolFromMolBlock(block, removeHs=False, sanitize=True)
        if molecule is None: raise ValueError(f'Cannot roundtrip residue {key}')
        if molecule.GetNumAtoms() != residue.GetNumAtoms(): raise ValueError('Atom count changed on roundtrip')
        for atom, original in zip(molecule.GetAtoms(), residue.GetAtoms()):
            atom.SetUnsignedProp('mapindex', original.GetUnsignedProp('mapindex'))
            # Molblocks omit residue identity. Preserve it for native water-aware
            # rules instead of accidentally turning HOH into UNK.
            if original.GetMonomerInfo(): atom.SetMonomerInfo(original.GetMonomerInfo())
        common = Residue(molecule)
        positions = common.GetConformer().GetPositions()
        atoms = [{'index': a.GetIdx(), 'parentIndex': a.GetUnsignedProp('mapindex'),
                  'atomicNumber': a.GetAtomicNum(), 'position': positions[a.GetIdx()].tolist()}
                 for a in common.GetAtoms()]
        result['components'][key] = {'molblock': block, 'atoms': atoms, 'residue': str(residue.resid),
                                    'graphSha256': hashlib.sha256(block.encode()).hexdigest()}
        cache[key] = common
        return common

    receptor = Chem.MolFromPDBFile(str(folder / 'receptor.pdb'), removeHs=True, sanitize=True)
    if receptor is None: raise ValueError('Native RDKit could not sanitize receptor.pdb')
    if receptor.GetNumAtoms() > 15000: raise ValueError('Resource exclusion: receptor exceeds 15000 atoms')
    ligands = []
    for path in sorted((folder / 'ligand_files').glob('*.sdf')):
        supplier = Chem.SDMolSupplier(str(path), removeHs=True, sanitize=True)
        for index, mol in enumerate(supplier):
            if mol is None:
                result['errors'].append({'stage': 'ligand-input', 'file': path.name, 'index': index, 'error': 'RDKit parse/sanitize failed'})
            else: ligands.append((f'{path.stem}:{index}', mol))
    if not ligands: result['errors'].append({'stage': 'selection', 'error': 'No usable ligand SDFs'})
    cache = {}
    for lane in result['lanes']:
        protein = Molecule.from_rdkit(Chem.AddHs(receptor, addCoords=True, addResidueInfo=True) if lane == 'explicit' else receptor)
        for ligand_id, mol in ligands:
            try:
                ligand_mol = Molecule.from_rdkit(Chem.AddHs(mol, addCoords=True) if lane == 'explicit' else mol)
                for li, lig in enumerate(ligand_mol):
                    ligand_key = f'{lane}:ligand:{ligand_id}:{li}'
                    common_lig = component(lig, ligand_key)
                    for ri, res in enumerate(protein):
                        dist2 = ((res.xyz[:, None, :] - lig.xyz[None, :, :]) ** 2).sum(axis=-1)
                        if float(dist2.min()) > 64: continue
                        protein_key = f'{lane}:protein:{ri}:{res.resid}'
                        try: common_prot = component(res, protein_key)
                        except Exception as exc:
                            result['errors'].append({'stage': 'serialization', 'component': protein_key, 'error': str(exc)})
                            continue
                        row = {'lane': lane, 'ligand': ligand_key, 'protein': protein_key, 'rules': {}}
                        for name, rule in rules.items():
                            try:
                                expected = [normalize(name, m) for m in rule.detect(common_lig, common_prot)]
                                row['rules'][name] = {'expected': expected}
                            except Exception as exc: row['rules'][name] = {'error': f'{type(exc).__name__}: {exc}'}
                            try:
                                row['rules'][name]['beforeSerialization'] = [normalize(name, m) for m in rule.detect(lig, res)]
                            except Exception as exc:
                                row['rules'][name]['beforeSerializationError'] = f'{type(exc).__name__}: {exc}'
                        result['pairs'].append(row)
            except Exception as exc:
                result['errors'].append({'stage': lane, 'ligand': ligand_id, 'error': f'{type(exc).__name__}: {exc}'})
    write(folder / 'native.json', result)
    print(folder.name, len(result['components']), 'components;', len(result['pairs']), 'pairs;', len(result['errors']), 'errors', flush=True)


def generate(manifest_path, out):
    if (out / 'manifest.json').exists(): raise FileExistsError('Choose a new output directory; do not overwrite a frozen corpus')
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(manifest_path.read_text())
    write(out / 'manifest.json', manifest)
    for case in manifest['cases']:
        folder = out / case['pdb']
        folder.mkdir()
        source = Path(case['source']).parent
        try:
            shutil.copyfile(source / 'receptor.pdb', folder / 'receptor.pdb')
            shutil.copytree(source / 'ligand_files', folder / 'ligand_files')
            files = [folder / 'receptor.pdb', *sorted((folder / 'ligand_files').glob('*.sdf'))]
            write(folder / 'inputs.json', {str(f.relative_to(folder)): sha(f) for f in files})
            env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
            with (folder / 'native.log').open('w') as log:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), '--case', str(folder)],
                               stdout=log, stderr=subprocess.STDOUT, env=env, check=True, timeout=180)
            print(case['pdb'], 'native complete', flush=True)
        except Exception as exc:
            write(folder / 'native.json', {'error': f'{type(exc).__name__}: {exc}'})
            print(case['pdb'], 'ERROR', str(exc), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--case', type=Path)
    args = parser.parse_args()
    if args.case: run_case(args.case)
    elif args.manifest and args.out: generate(args.manifest, args.out)
    else: parser.error('Supply --case or both --manifest and --out')
