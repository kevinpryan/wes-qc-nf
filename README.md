The goal of this pipeline is to take BAMs that have been produced by the nf-core/fastquorum pipeline and carry out QC to give us the following metrics:

- Chimeras
- Unmapped reads percent
- On target coverage
- Fold 80 base penalty
- False positive rate

To run the pipeline run the following:

```bash
nextflow run main.nf \
  --bams "/path/to/bamdir/*.bam" \
  --ref "/path/to/refdir/reference.fa" \
  --probes_bed "/path/to/bed-dir/IDT.Exomev2.probes.hg38.bed" \
  --targets_bed "/path/to/bed-dir/IDT.Exomev2.targets.hg38.bed" \
  --outdir "QC_Results" \
  -profile <singularity/cluster/.../institute>
```

Profiles taken from `nf-hlamajority`
