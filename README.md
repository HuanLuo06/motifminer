# MotifMiner

MotifMiner analyzes known ligand-binding positions in protein multiple sequence
alignments. It can analyze a user-provided MSA, use the original remote NCBI
workflow, or search any compatible local protein BLAST database.

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

The command writes `conservation.tsv` for inspection and `conservation.json`
for downstream software and the future web interface.

## Local BLAST workflow

`motifminer local` accepts a local database prefix produced by `makeblastdb`.
The database can contain GTDB, RefSeq, UniProt, or another protein collection;
MotifMiner does not hard-code database locations.

Domain/sensor query with canonical full-protein binding-site coordinates:

```bash
motifminer local --query-mode domain \
  --reference kind_full.fasta \
  --domain kind_sensor.fasta \
  --sites kind_sites.csv \
  --database /data/blast/gtdb_proteins \
  --database-name GTDB \
  --evalue 1e-5 --min-query-coverage 0.8 \
  --min-identity 0.25 --max-identity 0.95 \
  --threads 32 --output results/kind_gtdb
```

The domain sequence must occur exactly once in the full reference, and every
binding site must lie inside it. MotifMiner records the full-position to
domain-position to MSA-column conversion in `site_mapping_qc.tsv`.

Full-protein query, for example against RefSeq:

```bash
motifminer local --query-mode full \
  --reference dctb_full.fasta --sites dctb_sites.csv \
  --database /data/blast/refseq_protein --database-name RefSeq \
  --threads 24 --output results/dctb_refseq
```

UniProt is selected the same way by passing its local database prefix, for
example `--database /data/blast/uniprot_sprot --database-name UniProtKB-Swiss-Prot`.

Local BLAST emits aligned subject sequence data and coordinates directly. This
allows homologous regions to be extracted even when the database was built
without `makeblastdb -parse_seqids`. Outputs include raw and filtered BLAST
tables, homolog regions, alignment input, final alignment, mapping/QC,
conservation tables, and pipeline metadata.

## Legacy remote workflow

With NCBI BLAST+ and MAFFT installed, MotifMiner can submit a remote NCBI protein BLAST search, filter hits, retrieve complete protein sequences from NCBI,
align them, and calculate binding-site conservation:

```bash
motifminer run \
  --reference reference.fasta \
  --sites binding_sites.csv \
  --email you@example.org \
  --threads 8 \
  --output results
```

The stages are resumable. Existing BLAST, homolog FASTA, and alignment files in
the output directory are reused unless `--force` is supplied. NCBI remote
services are appropriate for occasional interactive analyses; a production web
server should use a maintained local database or a controlled job queue.

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
