# MotifMiner

MotifMiner analyzes known ligand-binding positions in a full-protein multiple
sequence alignment. Version 0.1 validates reference coordinates, maps them to
MSA columns, reports nongap amino-acid frequencies, and calculates exact,
automatic, extended, and user-defined grouped conservation.

## Current input

- A single-record reference protein FASTA.
- A binding-site CSV with `position`, `residue`, and optional
  `accepted_residues` columns.
- A full-protein aligned FASTA containing the unchanged reference sequence.

Positions are one-based reference-protein coordinates. Accepted residues may be
written as `YF`, `Y,F`, or `Y/F`.

## Run the demo

```bash
python -m motifminer learn \
  --reference examples/dctb_demo/reference.fasta \
  --sites examples/dctb_demo/binding_sites.csv \
  --msa examples/dctb_demo/alignment.fasta \
  --output demo_results
```

The command writes `conservation.tsv` for inspection and
`conservation.json` for downstream software and the future web interface.

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Interpretation

`nongap_n` counts sequences containing a standard amino acid at a binding-site
column. Amino-acid frequencies use `nongap_n` as their denominator. Gaps and
unknown characters are reported separately.

Automatic substitutions are deliberately conservative. Extended substitutions
are reported separately for user review because biochemical similarity does not
guarantee equivalent ligand binding.

