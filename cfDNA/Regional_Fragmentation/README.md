If this still appears as one paragraph, then **something above it is broken**.

---

## ✅ Correct Step 0 (SAFE VERSION)

Once the test works, paste **this EXACT block** (no extra indentation):

```md
## Step 0 — BAM preprocessing (coordinate sort + filtering)

### Goal
Generate coordinate-sorted BAM files containing high-confidence cfDNA fragments.

### Filtering criteria
- Proper pairs (`-f 2`)
- MAPQ ≥ 30 (`-q 30`)
- Remove secondary, QC-fail, duplicate, supplementary (`-F 3844`)

### Command
```bash
samtools view -bh -f 2 -q 30 -F 3844 sample.bam \
  | samtools sort -o sample.filter.coordsort.bam -

samtools index sample.filter.coordsort.bam
