Breseq reference sequences and run records
Prepared: 2026-09-20

Contents
- 149 study-run folders, each with:
  - data/reference.fasta
  - data/summary.json
  - output/summary.json
  - output/log.txt (the recorded breseq command line and reference path)
- Three GenBank references:
  - work/greencenter/s434226/plac/reference/data/all_ATEC_annotated_contigs_w_locustag.gbk
  - project/greencenter/Toprak_lab/shared/adam/adam_sequencing/reference/data/all_ATEC_annotated_contigs_w_locustag_w_genetags.gbk
  - project/greencenter/Toprak_lab/shared/adam/adam_sequencing/reference/data/sequence-4.gb

Availability and provenance audit
- All 149 run folders contained all four requested records (596 records total).
- All 149 output/log.txt files contain a breseq command line. The command records
  --polymorphism-prediction, --num-processors 16, the reference argument, and
  the original FASTQ input path.
- The 141 U00096 runs logged the now-absent source path
  /work/greencenter/s434226/plac/reference/data/sequence-4.gb. The included
  project copy is U00096.3 and has sequence content identical to those 141
  per-run reference.fasta files.
- The eight ATEC runs logged the now-absent source path
  /work/greencenter/s434226/plac/reference/data/all_ATEC_annotated_contigs_w_locustag_w_genetags.gbk.
  The included project copy has sequence content identical to their per-run
  reference.fasta files.
- No standalone shell, Slurm, PBS, job, or command script was found in the
  complete 40,754-file BioHPC out tree. The run logs contain no breseq version
  string; the summary JSON files retain the run options.

Integrity
- 599 source files in this archive were verified against BioHPC using SHA-256
  before the ZIP archive was made.
