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
    ImplicitHBAcceptor,
    ImplicitHBDonor,
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
from prolif.molecule import Molecule, sdf_supplier
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
    geometry: dict[str, float | list[float]] = {}
    for key, value in metadata.items():
        if key in common:
            continue
        if isinstance(value, (int, float, np.number)):
            geometry[key] = float(value)
        elif isinstance(value, (list, tuple, np.ndarray)) and all(
            isinstance(item, (int, float, np.number)) for item in value
        ):
            geometry[key] = [float(item) for item in value]
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


def residue_from_smiles(smiles: str, *, add_hydrogens: bool) -> Residue:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise RuntimeError(f"RDKit could not parse SMARTS parity molecule {smiles!r}")
    if add_hydrogens:
        molecule = Chem.AddHs(molecule)
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, atom in enumerate(molecule.GetAtoms()):
        conformer.SetAtomPosition(index, Point3D(float(index), 0.0, 0.0))
        atom.SetUnsignedProp("mapindex", index)
    molecule.AddConformer(conformer)
    return Residue(molecule)


def make_smarts_cases(
    pattern_key: str,
    queries: Chem.Mol | list[Chem.Mol],
    smiles_values: list[str],
    *,
    hydrogen_modes: tuple[bool, ...] = (False, True),
) -> list[dict[str, Any]]:
    query_list = queries if isinstance(queries, list) else [queries]
    output = []
    for smiles in smiles_values:
        for add_hydrogens in hydrogen_modes:
            component = residue_from_smiles(
                smiles,
                add_hydrogens=add_hydrogens,
            )
            expected_matches = [
                list(match)
                for query in query_list
                for match in component.GetSubstructMatches(query)
            ]
            output.append(
                {
                    "name": (
                        f"{pattern_key}:{smiles}:"
                        f"{'explicit-h' if add_hydrogens else 'implicit-h'}"
                    ),
                    "patternKey": pattern_key,
                    "pythonSmarts": [Chem.MolToSmarts(query) for query in query_list],
                    "component": serialize_component(component),
                    "expectedMatches": expected_matches,
                }
            )
    return output


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

    implicit_path = datapath / "implicitHbond"
    implicit_ligand = sdf_supplier(str(implicit_path / "1.D.sdf"))[0][0]
    implicit_protein_rdkit = Chem.MolFromPDBFile(
        str(implicit_path / "receptor.pdb"),
        removeHs=True,
    )
    if implicit_protein_rdkit is None:
        raise RuntimeError("RDKit could not parse the implicit-H receptor fixture")
    implicit_protein = Molecule.from_rdkit(implicit_protein_rdkit)
    implicit_tyr = implicit_protein["TYR167.B"]
    implicit_asp = implicit_protein["ASP95.A"]

    definitions: list[tuple[str, str, Interaction, Residue, Residue]] = [
        (
            "hydrophobic-edge",
            "Hydrophobic",
            Hydrophobic(),
            fixtures["benzene"],
            fixtures["edge"],
        ),
        (
            "hydrophobic-negative",
            "Hydrophobic",
            Hydrophobic(),
            fixtures["benzene"],
            fixtures["chlorine"],
        ),
        (
            "hydrophobic-inclusive-cutoff",
            "Hydrophobic",
            Hydrophobic(),
            bromine((0, 0, 0), 11),
            bromine((4.5, 0, 0), 29),
        ),
        (
            "hydrophobic-outside-cutoff",
            "Hydrophobic",
            Hydrophobic(),
            bromine((0, 0, 0), 13),
            bromine((4.5001, 0, 0), 31),
        ),
        (
            "hbond-acceptor",
            "HBAcceptor",
            HBAcceptor(),
            fixtures["acceptor"],
            fixtures["donor"],
        ),
        (
            "hbond-acceptor-negative",
            "HBAcceptor",
            HBAcceptor(),
            fixtures["acceptor_false"],
            fixtures["donor"],
        ),
        ("hbond-donor", "HBDonor", HBDonor(), fixtures["donor"], fixtures["acceptor"]),
        (
            "hbond-donor-negative",
            "HBDonor",
            HBDonor(),
            fixtures["donor"],
            fixtures["acceptor_false"],
        ),
        (
            "implicit-hbond-acceptor-tyr",
            "ImplicitHBAcceptor",
            ImplicitHBAcceptor(),
            implicit_tyr,
            implicit_ligand,
        ),
        (
            "implicit-hbond-acceptor-asp",
            "ImplicitHBAcceptor",
            ImplicitHBAcceptor(),
            implicit_asp,
            implicit_ligand,
        ),
        (
            "implicit-hbond-acceptor-negative",
            "ImplicitHBAcceptor",
            ImplicitHBAcceptor(),
            implicit_tyr,
            fixtures["chlorine"],
        ),
        (
            "implicit-hbond-donor",
            "ImplicitHBDonor",
            ImplicitHBDonor(),
            implicit_ligand,
            implicit_tyr,
        ),
        (
            "implicit-hbond-donor-negative",
            "ImplicitHBDonor",
            ImplicitHBDonor(),
            fixtures["chlorine"],
            implicit_tyr,
        ),
        (
            "halogen-acceptor",
            "XBAcceptor",
            XBAcceptor(),
            fixtures["xb_acceptor"],
            fixtures["xb_donor"],
        ),
        (
            "halogen-acceptor-xar-negative",
            "XBAcceptor",
            XBAcceptor(),
            fixtures["xb_acceptor_false_xar"],
            fixtures["xb_donor"],
        ),
        (
            "halogen-acceptor-axd-negative",
            "XBAcceptor",
            XBAcceptor(),
            fixtures["xb_acceptor_false_axd"],
            fixtures["xb_donor"],
        ),
        (
            "halogen-donor",
            "XBDonor",
            XBDonor(),
            fixtures["xb_donor"],
            fixtures["xb_acceptor"],
        ),
        (
            "halogen-donor-negative",
            "XBDonor",
            XBDonor(),
            fixtures["xb_donor"],
            fixtures["xb_acceptor_false_xar"],
        ),
        ("cationic", "Cationic", Cationic(), fixtures["cation"], fixtures["anion"]),
        (
            "cationic-negative",
            "Cationic",
            Cationic(),
            fixtures["cation"],
            fixtures["benzene"],
        ),
        ("anionic", "Anionic", Anionic(), fixtures["anion"], fixtures["cation"]),
        (
            "anionic-negative",
            "Anionic",
            Anionic(),
            fixtures["anion"],
            fixtures["benzene"],
        ),
        ("cation-pi", "CationPi", CationPi(), fixtures["cation"], fixtures["benzene"]),
        (
            "cation-pi-negative",
            "CationPi",
            CationPi(),
            fixtures["cation_false"],
            fixtures["benzene"],
        ),
        ("pi-cation", "PiCation", PiCation(), fixtures["benzene"], fixtures["cation"]),
        (
            "pi-cation-negative",
            "PiCation",
            PiCation(),
            fixtures["benzene"],
            fixtures["cation_false"],
        ),
        (
            "pi-stacking-face",
            "PiStacking",
            PiStacking(),
            fixtures["benzene"],
            fixtures["face"],
        ),
        (
            "pi-stacking-edge",
            "PiStacking",
            PiStacking(),
            fixtures["benzene"],
            fixtures["edge"],
        ),
        (
            "pi-stacking-negative",
            "PiStacking",
            PiStacking(),
            fixtures["benzene"],
            fixtures["chlorine"],
        ),
        (
            "metal-donor",
            "MetalDonor",
            MetalDonor(),
            fixtures["metal"],
            fixtures["chelator"],
        ),
        (
            "metal-donor-negative",
            "MetalDonor",
            MetalDonor(),
            fixtures["metal_false"],
            fixtures["chelator"],
        ),
        (
            "metal-acceptor",
            "MetalAcceptor",
            MetalAcceptor(),
            fixtures["chelator"],
            fixtures["metal"],
        ),
        (
            "metal-acceptor-negative",
            "MetalAcceptor",
            MetalAcceptor(),
            fixtures["chelator"],
            fixtures["metal_false"],
        ),
        (
            "vdw-contact",
            "VdWContact",
            VdWContact(),
            fixtures["benzene"],
            fixtures["edge"],
        ),
        (
            "vdw-negative",
            "VdWContact",
            VdWContact(),
            fixtures["acceptor"],
            fixtures["metal_false"],
        ),
    ]
    cases = [make_case(*definition) for definition in definitions]

    hydrophobic = Hydrophobic()
    explicit_hbond = HBAcceptor()
    implicit_hbond = ImplicitHBAcceptor()
    halogen = XBAcceptor()
    ionic = Cationic()
    cation_pi = CationPi()
    metal = MetalDonor()
    smarts_cases = [
        *make_smarts_cases(
            "Hydrophobic.lig_pattern",
            hydrophobic.lig_pattern,
            [
                "C",
                "C=[SH2]",
                "c1cscc1",
                "c1cocc1",
                "[*]SC",
                "[*]CC",
                "[*]C=C",
                "[*]=C=C",
                "[*]C(=C)C",
                "[*]C(C)C",
                "CS(C)(C)C",
                "FC(F)(F)F",
                "BrI",
                "C=O",
                "C=N",
                "CF",
                "Nc1ccccc1",
                "[*]C(C)(C)C",
            ],
        ),
        *make_smarts_cases(
            "HBAcceptor.lig_pattern",
            explicit_hbond.lig_pattern,
            [
                "O",
                "N",
                "[NH4+]",
                "C-N-C=O",
                "N-C=[SH2]",
                "[nH+]1ccccc1",
                "n1ccccc1",
                "n(C)1cccc1",
                "[nH]1cccc1",
                "[nH+]1nc[nH]c1",
                "c12c([nH]cc1)cccc2",
                "Nc1ccccc1",
                "C(=N)(N)N",
                "N#C",
                "o1cccc1",
                "[o+](C)1cccc1",
                "[*]-[N+](=O)-[O-]",
                "COC=O",
                "c1ccccc1Oc1ccccc1",
                "FC",
                "Fc1ccccc1",
                "FCF",
                "c1cc[n+](cc1)[O-]",
            ],
        ),
        *make_smarts_cases(
            "HBAcceptor.prot_pattern",
            explicit_hbond.prot_pattern,
            [
                "[OH2]",
                "[NH3]",
                "[NH4+]",
                "[SH2]",
                "O=C=O",
                "c1c[nH+]ccc1",
                "c1c[nH]cc1",
                "C-N-C=O",
            ],
            hydrogen_modes=(True,),
        ),
        *make_smarts_cases(
            "ImplicitHBAcceptor.lig_pattern",
            implicit_hbond.lig_pattern,
            ["O", "N", "[NH4+]", "n1ccccc1", "[nH]1cccc1", "C(=N)(N)N", "COC=O"],
        ),
        *make_smarts_cases(
            "ImplicitHBAcceptor.prot_pattern",
            implicit_hbond.prot_pattern,
            ["O", "N", "[NH4+]", "[SH2]", "c1c[nH+]ccc1", "c1c[nH]cc1", "C-N-C=O"],
        ),
        *make_smarts_cases(
            "XBAcceptor.lig_pattern",
            halogen.lig_pattern,
            ["[NH3]", "[NH+]C", "c1ccccc1", "C(=O)C", "C#N"],
        ),
        *make_smarts_cases(
            "XBAcceptor.prot_pattern",
            halogen.prot_pattern,
            ["CCl", "c1ccccc1Cl", "NCl", "c1cccc[n+]1Cl", "CC"],
        ),
        *make_smarts_cases(
            "Cationic.lig_pattern",
            ionic.lig_pattern,
            [
                "[NH4+]",
                "[Ca+2]",
                "CC(=[NH2+])N",
                "NC(=[NH2+])N",
                "c1cc[n+](cc1)[O-]",
                "C-N=[N+]=[N-]",
                "N#[N+]-[N-2]",
            ],
        ),
        *make_smarts_cases(
            "Cationic.prot_pattern",
            ionic.prot_pattern,
            [
                "[Cl-]",
                "CC(=O)[O-]",
                "CS(=O)[O-]",
                "CP(=O)[O-]",
                "c1cc[n+](cc1)[O-]",
                "C-N=[N+]=[N-]",
            ],
        ),
        *make_smarts_cases(
            "CationPi.cation",
            cation_pi.cation,
            [
                "[NH4+]",
                "[Ca+2]",
                "CC(=[NH2+])N",
                "NC(=[NH2+])N",
                "C-N=[N+]=[N-]",
                "N#[N+]-[N-2]",
            ],
        ),
        *make_smarts_cases(
            "CationPi.pi_ring",
            cation_pi.pi_ring,
            ["c1ccccc1", "c1cocc1", "C1CCCCC1", "CC"],
        ),
        *make_smarts_cases(
            "MetalDonor.lig_pattern",
            metal.lig_pattern,
            ["[Mg]", "[Zn]", "[Na+]", "O"],
        ),
        *make_smarts_cases(
            "MetalDonor.prot_pattern",
            metal.prot_pattern,
            [
                "O",
                "N",
                "[NH+]",
                "N-C=[SH2]",
                "[nH+]1ccccc1",
                "Nc1ccccc1",
                "o1cccc1",
                "COC=O",
                "[*]-[N+](=O)-[O-]",
            ],
        ),
    ]
    return {
        # Editable and CI builds can report 0.0.0 for this same source tree.
        # Keep oracle provenance tied to the rule version instead of packaging.
        "generatedBy": "ProLIF 2.2.1",
        "pythonRdkitVersion": rdkit.__version__,
        "rules": sorted({case["rule"] for case in cases}),
        "cases": cases,
        "smartsCases": smarts_cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(generate(), indent=2) + "\n")


if __name__ == "__main__":
    main()
