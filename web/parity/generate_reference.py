"""Generate browser-port parity fixtures with upstream Python ProLIF.

This is a development oracle, not browser runtime code. Run it from the ProLIF
repository so the checked-out Python implementation is the source of truth.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import prolif
import rdkit
from MDAnalysis import Universe
from MDAnalysis.topology.guessers import guess_atom_element
from rdkit import Chem
from rdkit.Geometry import Point3D

from prolif.datafiles import datapath
from prolif.interactions import Hydrophobic
from prolif.interactions.utils import get_mapindex
from prolif.molecule import Molecule
from prolif.residue import Residue


def from_mol2(filename: str) -> Residue:
    """Load the same interaction fixtures used by ProLIF's Python tests."""
    universe = Universe(str(datapath / filename))
    elements = [guess_atom_element(name) for name in universe.atoms.names]
    universe.add_TopologyAttr("elements", np.asarray(elements, dtype=object))
    universe.atoms.types = np.asarray(
        [atom_type.upper() for atom_type in universe.atoms.types],
        dtype=object,
    )
    return Molecule.from_mda(universe, force=True)[0]


def bromine(position: tuple[float, float, float], parent_index: int) -> Residue:
    molecule = Chem.MolFromSmiles("[Br]")
    if molecule is None:  # pragma: no cover - defensive RDKit boundary
        raise RuntimeError("RDKit could not construct bromine fixture")
    conformer = Chem.Conformer(1)
    conformer.SetAtomPosition(0, Point3D(*position))
    molecule.AddConformer(conformer)
    molecule.GetAtomWithIdx(0).SetUnsignedProp("mapindex", parent_index)
    return Residue(molecule)


def serialize_component(residue: Residue) -> dict[str, Any]:
    atoms = []
    conformer = residue.GetConformer()
    for atom in residue.GetAtoms():
        index = atom.GetIdx()
        point = conformer.GetAtomPosition(index)
        atoms.append(
            {
                "index": index,
                "parentIndex": get_mapindex(residue, index),
                "atomicNumber": atom.GetAtomicNum(),
                "position": [point.x, point.y, point.z],
            }
        )
    return {
        "molblock": Chem.MolToMolBlock(residue),
        "atoms": atoms,
    }


def normalize_metadata(metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        "interaction": "Hydrophobic",
        "indices": {
            "ligand": list(metadata["indices"]["ligand"]),
            "protein": list(metadata["indices"]["protein"]),
        },
        "parentIndices": {
            "ligand": list(metadata["parent_indices"]["ligand"]),
            "protein": list(metadata["parent_indices"]["protein"]),
        },
        "distance": metadata["distance"],
    }


def make_case(name: str, ligand: Residue, protein: Residue) -> dict[str, Any]:
    expected = [
        normalize_metadata(metadata)
        for metadata in Hydrophobic().detect(ligand, protein)
    ]
    return {
        "name": name,
        "ligand": serialize_component(ligand),
        "protein": serialize_component(protein),
        "expected": expected,
    }


def generate() -> dict[str, Any]:
    benzene = from_mol2("benzene.mol2")
    cases = [
        make_case("benzene-edge-to-face", benzene, from_mol2("edgetoface.mol2")),
        make_case("benzene-face-to-face", benzene, from_mol2("facetoface.mol2")),
        make_case("benzene-chlorine-negative", benzene, from_mol2("chlorine.mol2")),
        make_case("benzene-bromine", benzene, from_mol2("bromine.mol2")),
        make_case("benzene-anion-negative", benzene, from_mol2("anion.mol2")),
        make_case("benzene-cation-negative", benzene, from_mol2("cation.mol2")),
        make_case("inclusive-cutoff", bromine((0, 0, 0), 11), bromine((4.5, 0, 0), 29)),
        make_case(
            "outside-cutoff",
            bromine((0, 0, 0), 13),
            bromine((4.5001, 0, 0), 31),
        ),
    ]
    return {
        "generatedBy": f"ProLIF {prolif.__version__}",
        "pythonRdkitVersion": rdkit.__version__,
        "rule": "Hydrophobic",
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(generate(), indent=2) + "\n")


if __name__ == "__main__":
    main()
