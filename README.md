# DCM comorbidity PRS
Cardiometabolic polygenic risk in familial dilated cardiomyopathy

This repository contains the code for a study of polygenic risk scores (PRS) for cardiometabolic co-morbidities in
familial dilated cardiomyopathy (DCM). Briefly, odds ratios per standard deviation of each PRS were estimated by
logistic regression in a discovery cohort (VCCRI) and a validation cohort (AGHA) of European-ancestry
DCM probands, each compared with a separate half of the Medical Genome Reference Bank (MGRB) controls.
All probands were also compared with 1000 Genomes non-Finnish Europeans as a sensitivity analysis, and
discrimination was assessed by the area under the curve (AUC) in all probands and all MGRB controls.

This study is described in: Thompson JM, Johnson R, Young PE, Soka MJ, Ohanian M, Fatkin D, Giannoulatou E.
Cardiometabolic polygenic risk in familial dilated cardiomyopathy. Under review.

Genomic data for this study were accessed through Australian Genomics, funded by the NHMRC (1113531, 2000001) and the
MRFF, administered by the Murdoch Children's Research Institute. The Australian Genomics Cardiovascular Disorders
Flagship was funded by the MRFF Genomic Health Futures Mission (EPCD000028).

## Scoring

PRS were calculated with PLINK 2 (v2.0.0) from published score files in GRCh38 (PGS Catalog harmonised files,
and the HUNT score for cardiorespiratory fitness). For each score, cases and
controls were restricted to the variants present in both datasets and in the score file, and the same variant QC was
applied to each:

plink2 --pfile cohort --extract shared_variants.txt --rm-dup exclude-mismatch --maf 0.01 --geno 0.05 --mind 0.05 --make-pgen --out cohort_shared

Both datasets were then scored on the variants that passed QC in both:

plink2 --pfile cohort_shared --extract scoring_ids.txt --score score_shared.txt 1 2 3 cols=+scoresums --out prs_cohort_score

## Analysis

The DCM_comorbidity_prs_analysis.R script was used to produce the figures and statistics in the research letter from the PLINK 2 output. The script requires R with ggplot2, patchwork and pROC (analysed with R 4.5.2, ggplot2 4.0.2, patchwork 1.3.2, pROC 1.19.0.1).

Genotype and score data are not publicly available due to patient privacy and consent requirements.
