import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
} from "../types.js";
import { invertedMetadata } from "./metadata.js";

export class InvertedRule implements InteractionRule {
  constructor(
    readonly name: string,
    private readonly source: InteractionRule,
  ) {}

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    return this.source
      .detect(protein, ligand)
      .map((metadata) => invertedMetadata(this.name, metadata));
  }
}
