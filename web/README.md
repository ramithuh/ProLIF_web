# ProLIF Web

`@ramithuh/prolif-web` is an **unofficial**, browser-native TypeScript port of
the deterministic interaction rules in [ProLIF](https://github.com/chemosim-lab/ProLIF).
It is intended to let applications such as Weaver calculate interaction
fingerprints locally, without a Python API or server-side chemistry worker.

This package is an interaction engine, not a viewer. PocketScope and Weaver
remain separate applications that consume its results.

## Compatibility target

The initial target is ProLIF 2.2.x. Upstream Python ProLIF remains the reference
implementation and parity oracle. Rules are ported individually, with fixtures
that compare atom participation and geometry against upstream behavior.

The first milestone provides:

- stable molecule, feature-match, and interaction metadata contracts;
- browser-safe vector geometry;
- the generic distance interaction primitive used by rules such as
  `Hydrophobic` and `Cationic`;
- an RDKit.js SMARTS adapter.

Higher-order angle, ring, and van der Waals rules will be added rule by rule,
with explicit compatibility metadata rather than silently claiming full
ProLIF parity.

## Development

```bash
cd web
npm install
npm test
npm run typecheck
npm run build
```

## Status and naming

This project is not affiliated with or endorsed by the ProLIF or RDKit teams.
“ProLIF” is used only to describe the upstream project and compatibility goal.
See [THIRD_PARTY_NOTICES.md](./THIRD_PARTY_NOTICES.md) for provenance and
license details.
