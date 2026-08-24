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
from prolif.interactions import (
    Anionic,
    Cationic,
    CationPi,
    HBAcceptor,
    HBDonor,
    Hydrophobic,
    MetalAcceptor,
    MetalDonor,
    PiCation,
    PiStacking,
    VdWContact,
    XBAcceptor,
    XBDonor,
)
from prolif.interactions.base import Interaction
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
    return {"molblock": Chem.MolToMolBlock(residue), "atoms": atoms}


def normalize_metadata(
    interaction_name: str,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    common = {"indices", "parent_indices", "distance"}
    geometry = {
        key: float(value)
        for key, value in metadata.items()
        if key not in common and isinstance(value, (int, float, np.number))
    }
    normalized: dict[str, Any] = {
        "interaction": interaction_name,
        "indices": {
            "ligand": list(metadata["indices"]["ligand"]),
            "protein": list(metadata["indices"]["protein"]),
        },
        "parentIndices": {
            "ligand": list(metadata["parent_indices"]["ligand"]),
            "protein": list(metadata["parent_indices"]["protein"]),
        },
        "distance": float(metadata["distance"]),
    }
    if geometry:
        normalized["geometry"] = geometry
    return normalized


def make_case(
    name: str,
    interaction_name: str,
    rule: Interaction,
    ligand: Residue,
    protein: Residue,
) -> dict[str, Any]:
    expected = [
        normalize_metadata(interaction_name, metadata)
        for metadata in rule.detect(ligand, protein)
    ]
    return {
        "name": name,
        "rule": interaction_name,
        "ligand": serialize_component(ligand),
        "protein": serialize_component(protein),
        "expected": expected,
    }


def generate() -> dict[str, Any]:
    fixtures = {
        name: from_mol2(filename)
        for name, filename in {
            "benzene": "benzene.mol2",
            "edge": "edgetoface.mol2",
            "face": "facetoface.mol2",
            "chlorine": "chlorine.mol2",
            "bromine": "bromine.mol2",
            "anion": "anion.mol2",
            "cation": "cation.mol2",
            "cation_false": "cation_false.mol2",
            "acceptor": "acceptor.mol2",
            "acceptor_false": "acceptor_false.mol2",
            "donor": "donor.mol2",
            "xb_acceptor": "xbond_acceptor.mol2",
            "xb_acceptor_false_xar": "xbond_acceptor_false_xar.mol2",
            "xb_acceptor_false_axd": "xbond_acceptor_false_axd.mol2",
            "xb_donor": "xbond_donor.mol2",
            "metal": "metal.mol2",
            "metal_false": "metal_false.mol2",
            "chelator": "ligand.mol2",
        }.items()
    }

    definitions: list[tuple[str, str, Interaction, Residue, Residue]] = [
        ("hydrophobic-edge", "Hydrophobic", Hydrophobic(), fixtures["benzene"], fixtures["edge"]),
        ("hydrophobic-negative", "Hydrophobic", Hydrophobic(), fixtures["benzene"], fixtures["chlorine"]),
        ("hydrophobic-inclusive-cutoff", "Hydrophobic", Hydrophobic(), bromine((0, 0, 0), 11), bromine((4.5, 0, 0), 29)),
        ("hydrophobic-outside-cutoff", "Hydrophobic", Hydrophobic(), bromine((0, 0, 0), 13), bromine((4.5001, 0, 0), 31)),
        ("hbond-acceptor", "HBAcceptor", HBAcceptor(), fixtures["acceptor"], fixtures["donor"]),
        ("hbond-acceptor-negative", "HBAcceptor", HBAcceptor(), fixtures["acceptor_false"], fixtures["donor"]),
        ("hbond-donor", "HBDonor", HBDonor(), fixtures["donor"], fixtures["acceptor"]),
        ("hbond-donor-negative", "HBDonor", HBDonor(), fixtures["donor"], fixtures["acceptor_false"]),
        ("halogen-acceptor", "XBAcceptor", XBAcceptor(), fixtures["xb_acceptor"], fixtures["xb_donor"]),
        ("halogen-acceptor-xar-negative", "XBAcceptor", XBAcceptor(), fixtures["xb_acceptor_false_xar"], fixtures["xb_donor"]),
        ("halogen-acceptor-axd-negative", "XBAcceptor", XBAcceptor(), fixtures["xb_acceptor_false_axd"], fixtures["xb_donor"]),
        ("halogen-donor", "XBDonor", XBDonor(), fixtures["xb_donor"], fixtures["xb_acceptor"]),
        ("cationic", "Cationic", Cationic(), fixtures["cation"], fixtures["anion"]),
        ("cationic-negative", "Cationic", Cationic(), fixtures["cation_false"], fixtures["anion"]),
        ("anionic", "Anionic", Anionic(), fixtures["anion"], fixtures["cation"]),
        ("cation-pi", "CationPi", CationPi(), fixtures["cation"], fixtures["benzene"]),
        ("cation-pi-negative", "CationPi", CationPi(), fixtures["cation_false"], fixtures["benzene"]),
        ("pi-cation", "PiCation", PiCation(), fixtures["benzene"], fixtures["cation"]),
        ("pi-stacking-face", "PiStacking", PiStacking(), fixtures["benzene"], fixtures["face"]),
        ("pi-stacking-edge", "PiStacking", PiStacking(), fixtures["benzene"], fixtures["edge"]),
        ("metal-donor", "MetalDonor", MetalDonor(), fixtures["metal"], fixtures["chelator"]),
        ("metal-donor-negative", "MetalDonor", MetalDonor(), fixtures["metal_false"], fixtures["chelator"]),
        ("metal-acceptor", "MetalAcceptor", MetalAcceptor(), fixtures["chelator"], fixtures["metal"]),
        ("vdw-contact", "VdWContact", VdWContact(), fixtures["benzene"], fixtures["edge"]),
        ("vdw-negative", "VdWContact", VdWContact(), fixtures["acceptor"], fixtures["metal_false"]),
    ]
    cases = [make_case(*definition) for definition in definitions]
    return {
        "generatedBy": f"ProLIF {prolif.__version__}",
        "pythonRdkitVersion": rdkit.__version__,
        "rules": sorted({case["rule"] for case in cases}),
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
