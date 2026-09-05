# Real-structure native ProLIF differential tests

These are **strict, opt-in conformance tests**, separate from the ordinary
passing unit suite. Native-unsupported cases still make the full audit exit nonzero. They
are not `test.fails` tests and do not bless the current mismatches.

## Fast frozen regressions

From `web/`, with existing npm dependencies:

```bash
npx vitest run --config parity/plinder/regressions.config.ts
```

Three small fixtures from the random local sample reproduce an implicit-H-bond
neighbor-angle ordering mismatch (`7GMZ`), plane-angle mismatch (`2WYI`), and
named-water handling (`1WNU`).
They require no Python, local Plinder installation or network. All three now pass.
See [post-fix findings](FIXES-20260904.md) for the supported-check results and limits.

## Full corpus

```bash
/home/ruh/miniconda3/envs/xymol/bin/python parity/plinder/generate.py \
  --manifest /home/ruh/research/PhD/weaver/parity-runs/plinder-20260904/manifest.json \
  --out /absolute/path/to/new-corpus
PROLIF_PARITY_DIR=/absolute/path/to/new-corpus \
  npx vitest run --config parity/plinder/vitest.config.ts
```

The manifest contains 16 seeded-random distinct PDB entries, one available
system per entry, drawn from locally unpacked Plinder. The generator reads all
ligand SDFs and `receptor.pdb` from those systems. It does not download data,
filter to easy ligands or silently replace failures. Each native case runs
sequentially with a 180-second timeout and numerical-library threads limited
to one. Receptors above 15,000 atoms are retained as resource exclusions.
Existing corpus manifests cannot be overwritten.

Python ProLIF 2.2.1 is the reference. Each ligand/protein residue pair within
8 Angstrom of any atom is tested under all 15 ported rules, plus the CSD-radius
VdWContact configuration used by Weaver. This spatial selection is shared by
both engines; it is not a test of high-level residue selection. Rules with no
positive hits are visible in the coverage table, not certified by zeroes.

Two molecular input lanes are generated: heavy/implicit H, and explicit H added
once by native RDKit. Both Python and Web consume the same serialized molecular
graphs, coordinates, atom mappings and H states. The explicit lane does not test
independently generated hydrogen coordinates. Native residue identity is retained
for water-aware rules; the Web component API now accepts that context as `residueName`.

The audit compares event **multisets**, exact ordered local and parent atom
indices, and all reported geometry (absolute tolerance 1e-5). Neighbor-angle
array permutations are reported as metadata differences, not missing events.
Empty and omitted geometry objects are equivalent. Missing fields, extra fields,
nonfinite numbers and thrown errors cannot become matches.

Native errors are also exercised in Web. Both-sided errors are reported as
unsupported checks, excluded from successful comparison counts and still fail
the full-coverage gate. One-sided errors are behavioral failures. Native results
before residue serialization are compared separately (geometry tolerance 0.01)
to expose graph/context changes rather than attributing them to the port.

`PROLIF_PARITY_REPORT_ONLY=1` permits diagnostic report generation without a
failing process exit; it does **not** mean the parity test passed.

The saved audit records versions, source/input/WASM hashes, source snapshots,
all native expected events, per-rule coverage and detailed failures. Python
RDKit is 2025.03.6 and Web RDKit is 2025.03.4 in this run. They are not identical
versions; version-sensitive discrepancies require follow-up.

## Boundaries

This tests the implemented **rule engine and RDKit.js adapter**, not Weaver's
CCD download/preparation, ligand selection, UI, or ProLIF's trajectory/high-level
fingerprint machinery. Plain PDB-derived receptors have limited charge assignment;
this random sample has no positive salt-bridge/cation-pi/halogen-bond coverage.
Use the existing synthetic/native-fixture suite as complementary rule coverage,
not as a substitute for this real-structure audit or future charged-structure tests.

The first results are in [FINDINGS-20260904.md](FINDINGS-20260904.md).
