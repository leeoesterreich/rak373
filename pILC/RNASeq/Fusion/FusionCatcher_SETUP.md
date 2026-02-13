*FusionCatcher (v1.33, Singularity) – Troubleshooting & Fix Log*
```txt
Environment

FusionCatcher version: 1.33

Container:
/software/rhel9/manual/install/fusioncatcher/fusioncatcher-1.33.sif

HPC: Slurm-based cluster

Input: paired-end FASTQ

Database: human_v102

DB path:
/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/FusionCatcher_db/current
```
❗ Issue 1: faToTwoBit OpenSSL Error
🔴 Error Observed

During runtime:

```txt
faToTwoBit: error while loading shared libraries:
libssl.so.3: cannot open shared object file: No such file or directory
```

Even testing inside container:
```bash
singularity exec fusioncatcher-1.33.sif faToTwoBit
```
Produced:
```txt
error while loading shared libraries: libssl.so.3
```

# Root Cause

The faToTwoBit binary shipped inside the FusionCatcher SIF was dynamically linked against OpenSSL 3.

The container did not include libssl.so.3 / libcrypto.so.3.

Attempting to bind host OpenSSL libraries resulted in:

segmentation faults

ABI instability

Therefore, the container's faToTwoBit binary was broken for this system.

❗ Issue 2: Segmentation Fault After Binding OpenSSL

After binding:
```bash
--bind /lib64/libssl.so.3:/lib64/libssl.so.3
--bind /lib64/libcrypto.so.3:/lib64/libcrypto.so.3
```
faToTwoBit produced:
```java
Segmentation fault (core dumped)
```
Root Cause

Mixing host OpenSSL libraries with container glibc/loader caused ABI mismatch.

Injecting host system libraries into container is unsafe unless built against identical runtime.

*Final Stable Solution*
✔ Replace broken faToTwoBit with UCSC official binary

*Downloaded working UCSC binary:*
```bash
mkdir -p /ix1/alee/LO_LAB/Personal/Rahul/tools/ucsc
cd /ix1/alee/LO_LAB/Personal/Rahul/tools/ucsc

wget -O faToTwoBit http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x faToTwoBit

./faToTwoBit # Check
```
Works without OpenSSL dependency.

Override binary inside container at runtime
In Slurm script:
```bash
UCSC_FATOTWOBIT="/ix1/alee/LO_LAB/Personal/Rahul/tools/ucsc/faToTwoBit"
```
Run FusionCatcher with override:
```bash
singularity exec --cleanenv \
  --bind /ix1:/ix1 \
  --bind /tmp:/tmp \
  --bind "$UCSC_FATOTWOBIT":/usr/local/bin/faToTwoBit \
  "$SIF" \
  fusioncatcher.py \
    -d "$DB" \
    -i "$SAMPLE_INDIR" \
    -o "$SAMPLE_OUTDIR" \
    -T "$TMPDIR" \
    -p "$THREADS" \
    -Z
```


