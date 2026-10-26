# PDO organoid RNA-seq: fusion detection (nf-core/rnafusion)

Gene-fusion detection on bulk RNA-seq of patient-derived organoids (PDO), Lee-Oesterreich Lab.
Run on the Pitt CRC HTC cluster (Slurm, `--clusters=htc`) with nf-core/rnafusion 4.1.3.
By: Rahul Kumar.

## Pipeline in short
- Pipeline: nf-core/rnafusion **4.1.3**, Nextflow **25.10.4**, Singularity containers
- Callers (`--tools all`): Arriba, STAR-Fusion, FusionCatcher; FusionInspector checks read support;
  fusion-report merges the calls
- `--no_cosmic`: COSMIC is not used, so no login or password is needed anywhere
- Reference: human GRCh38, GENCODE 46, in `/ix1/alee/LO_LAB/General/References/rnafusion_4.1.3`
  (record of how it was downloaded: `Script/download_rnafusion_refs.sh`)
- Data: paired-end RNA-seq; all libraries so far are reverse-stranded (checked with Salmon)

## Batches
| Batch | Samples | Results |
|---|---|---|
| 1 (20220310_IPMBO_Jade_Organoid_RNAseq) | 55 | `20220310_IPMBO_Jade_Organoid_RNAseq/results/` |
| 2 (20240315_Organoid_RNAseq) | 41 | `20240315_Organoid_RNAseq/results/` |
| 3 (20260921_OrganoidBulkRNABatch3_Alex_Brent_Daniel) | 16 | `20260921_OrganoidBulkRNABatch3/results/` |

`samplesheet_fusion.csv` always holds the batch currently being run. When the batch is finished,
the sheet is moved into that batch's folder.

## Folder layout
- `Script/` final scripts (see below)
- `salmon_index_GRCh38_114/` Salmon index, built once and reused by every batch
- `logs/` all Slurm logs
- `20220310_IPMBO_Jade_Organoid_RNAseq/` key results of batch 1
- `20240315_Organoid_RNAseq/` key results of batch 2 (sample names start with a number, which the
  pipeline rejects, so the prefix `CHP_` was added)
- `20260921_OrganoidBulkRNABatch3/` key results of batch 3

Large pipeline files (work/, results) live on **/vast** (`/vast/alee/Rahul/PDO/RNAFusion_4.1.3`)
because /ix1 is almost full. Only the key results are copied to /ix1, then /vast is cleared.

## Scripts (in `Script/`)
| Script | Purpose |
|---|---|
| `build_salmon_index.sh` | Builds the Salmon index (once per species) |
| `check_strand_salmon.sh` | Salmon `-l A` on 1M read pairs per sample; writes `strand_summary.tsv` for the batch |
| `summarise_strand.sh` | Builds `strand_summary.tsv` (called automatically by the script above) |
| `make_samplesheet.sh` | Builds `samplesheet_fusion.csv` for one batch (see the notes below) |
| `download_rnafusion_refs.sh` | Record of how the rnafusion references were downloaded |
| `run_rnafusion_v4.sh` launcher, `rnafusion_v4.config` Nextflow settings |
| `copy_results_to_ix1.sh` | Copies the key results from /vast to /ix1 (no BAM, CRAM or FASTQ files); set `DEST` for each batch |
| `clear_vast.sh` | Deletes work/ and results/ on /vast after a verified copy; keeps the Singularity caches |

## How to run a new batch (step by step)
All commands from this folder.

1. **Strandedness.** Set the FASTQ folder and output name at the top of `Script/check_strand_salmon.sh`
   (`FASTQ_DIR`, `OUT_NAME`), then `sbatch Script/check_strand_salmon.sh`.
   Read `<OUT_NAME>/strand_summary.tsv`: `reverse` or `forward` if the fraction is above 0.8.
   Use a different Salmon index (`CDNA_GZ`, `INDEX`) if the samples are not human.
2. **Samplesheet.** Set the three defaults at the top of `make_samplesheet.sh`, then `bash make_samplesheet.sh`.
   It takes the strandedness from the summary, never overwrites an existing sheet, and stops
   on duplicate names. Columns: `sample,fastq_1,fastq_2,strandedness`.
   If sample names are only digits (108, 120 ...), run `PREFIX=CHP_ bash make_samplesheet.sh`:
   nf-core/rnafusion reads digit-only names as integers and rejects them, because sample names must be text.
3. **Run the pipeline.** Set the `SHEET` and `OUTDIR` defaults in `run_rnafusion_v4.sh`, then
   `sbatch run_rnafusion_v4.sh`. Check with `squeue --clusters=htc -u $USER`.
   The head job has a 72 h limit. If it ends early, just resubmit: `-resume` continues.
4. **Copy to /ix1.** Set `SRC` and `DEST` in `copy_results_to_ix1.sh` and run it. Check that every
   folder was copied before deleting anything.
5. **Clear /vast.** Run `clear_vast_batch3.sh` (set `DEST` for the batch) only after the copy is checked.
   Check the space with `crc-quota`.
6. **Move the samplesheet** into the batch folder.
7. **Consensus.** After all batches are done, run `build_consensus_v4.py` on all copied results
  .

## Consensus fusion analysis (all batches together)
Script: `Script/build_consensus_v4.py`. Results: `consensus_fusion/`.
Arriba, STAR-Fusion and FusionCatcher each have their own false positives, so their calls are merged
across all batches to find the fusions worth following up.

Approach:
1. Read every caller’s output into one format (genes, breakpoints, supporting reads, tags).
2. Calls with both breakpoints within 10 bp are the same event; n_tools counts how many callers agree.
3. FusionInspector read support is added as extra evidence (it re-checks the reads, so it is not a fourth vote).
4. Likely artifacts are flagged, not removed: readthrough or adjacent genes (same chromosome, closer than 100 kb), paralog and duplicated gene families, FusionCatcher normal tags (healthy, banned, gtex), unannotated genes.
5. Known fusions are labelled in known_fusion_db with FusionCatcher’s database tags (known, tcga, 18cancers, mitelman, ccle3, tcga3). This is a label only.
6. Shortlist = events called by at least 2 callers (--min_tools). Flags, tags, FusionInspector reads and the number of samples are kept as columns for the reader to judge..
Run ONCE for all batches together (list the results folder of each batch after --results; sample names must be unique across batches):

    sbatch --clusters=htc --cpus-per-task=8 --mem=16G --time=02:00:00 --job-name=consensus \
      --output=logs/consensus_%j.out --error=logs/consensus_%j.err \
      --wrap='module load anaconda3/2025.7.0-2-python_3.11
              cd /ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion
              python Script/build_consensus_v4.py \
                --results 20220310_IPMBO_Jade_Organoid_RNAseq/results 20240315_Organoid_RNAseq/results 20260921_OrganoidBulkRNABatch3/results \
                --outdir consensus_fusion'

For a new batch, add its results folder to `--results`.

Output files in `consensus_fusion/`:
- `fusions_all_calls.tsv` every raw call, one row per call per tool
- `fusions_consensus_per_sample.tsv` merged events per sample, with tools, support and flags
- `fusions_event_recurrence.tsv` each event across samples and batches
- `fusions_shortlist.tsv` the candidates to look at first

How the files are made (each one is built from the one above it):

```
all_calls            every call, every tool
   |  merge calls within 10 bp, add FusionInspector and flags
   v
consensus_per_sample    one event per sample
   |-- group the same event over all samples --> event_recurrence
   |-- keep rows passing the 4 rules ----------> shortlist
```

To check whether a fusion was ever detected, search `fusions_all_calls.tsv` (nothing is removed there).
A fusion missing from the shortlist can still be present in the other three files.

A shortlisted fusion is a candidate, not a proof. Next steps are checking the reads in IGV (FusionInspector
BAMs), the reading frame, and RT-PCR validation. Check whether samples sharing a fusion come from the same patient.
