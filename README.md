The goal of this pipeline is to take BAMs that have been produced by the nf-core/fastquorum pipeline and carry out QC to give us the following metrics:

- Chimeras
- Unmapped reads percent
- On target coverage
- Fold 80 base penalty
- False positive rate

To run the pipeline run the following:

```bash
nextflow run main.nf \
  --bams "s3://healed-data/pilot-run/wes-qc/bams/*.bam" \
  --ref "s3://healed-data/lens-v1.7-dev-healed/references/cloud/Homo_sapiens.assembly38.fa" \
  --probes_bed "s3://healed-data/pilot-run/wes-qc/references/IDT.Exomev2.probes.hg38.bed" \
  --targets_bed "s3://healed-data/pilot-run/wes-qc/references/IDT.Exomev2.targets.hg38.bed" \
  --outdir "QC_Results"
```

Profiles taken from `nf-hlamajority`
