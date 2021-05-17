# How to run on JLIM 2.0 sample data
JLIM 2.0 provides python scripts to generate permutation data for meta-analyzed secondary trait cohort. This allows to directly account for the subcohort-level variation in data processing and QC in meta-analyzed secondary trait data.

## Example data
We provide example data hosted in the Sunyaev lab website: [Example Data](http://genetics.bwh.harvard.edu/wiki/sunyaevlab/_media/jlim_2_example_revised.tar.gz). Please untar the example data in the same folder that the code was untared. The `examples_2.0/` folder should be in parallel with the `bin/` and `R/` directories. This dataset contains the following simulated data:
* list of index SNPs for height GWAS (Height/Height_indexSNP.tsv) 
* Height GWAS summary statistics 
* peaks file containing the midpoints of two sample genomic intervals to analyze
* reference LD files for two sample intervals

In addition, the example dataset includes the following files for each of three simulated cohorts (A, B, and C): 
* bimbam genotype file
* map file 
* phenotype file
* sample file 
* covariate file 

### Bimbam file
One genotype file in BIMBAM format is needed per cohort, per chromosome analyzed, for the secondary trait. These files should have no header and should be gzipped. The separator used in the file should be provided (see bimbam.separator). *Note: Currently, JLIM does not allow missing genotypes in this file.* 

### Map file 
The map (or info) file is a chromosome- and cohort-specific file with information on all variants present in the corresponding bimbam file. It should have no header and it should be gzipped. The number of lines in the map file and bimbam file should match exactly and be in the same order, with each line corresponding to a variant. The separator used in the file should be provided (see map.separator).

### Peaks file 
This file contains a comma-separated list of genomic positions to be analyzed in the corresponding chromosome. These will be at the center of intervals of size (in bp) “interval.int”, which will be analyzed using JLIM.

### Samples file 
This is a file with no header with the names of all samples in the bimbam file. This file has only one column, and there should be one sample in each row. The order of samples should match that of bimbam file. The sample names will be used to find the appropriate phenotypes and covariates. 

### Secondary phenotype file
This file has sample names in the first column and their corresponding phenotypes in the second column. Phenotypes that correspond to data sets where individual level data is provided are called secondary phenotypes, as opposed to primary phenotypes, where only summary statistics are needed. They should either be quantitative values or “nan”, which will exclude the sample from the analysis. The file should have a header and should be tab separated.

### Covariate file
This is a tab separated file which has a header with covariate names. Rows are samples, columns are covariates. The first column (column 0) contains sample names. There should be no missing values (‘nan’) in any of the covariates, otherwise samples will be excluded. You are free to include samples and covariates not used in the regression in this file, as the samples used will be specified in the phenotype file, and the covariates used will be specified in the argument covariates.list.

### indexSNP file 
The indexSNP file is a tab-delimiated file with five columns:
- CHR: chromosome
- SNP: SNP ID of index SNP
- BP: base pair position of index SNP
- STARTBP: start of interval around index SNP (bp)
- ENDBP: end of interval around index SNP (bp)
  
*CHR*.*STARTBP*.*ENDBP* combination will be used as an interval identifier to locate files of primary association statistics, secondary association statistics, and reference LD. In each interval, the most associated SNP will be automatically picked based on primary association p-values, and then the analysis window will be set up to +/- 100kb around the most associated SNP. The JLIM analysis window will not extend over the original interval specificied in the indexSNP file. The exact bp position of index SNP does not matter for JLIM as it will pick the most associated SNP within the specified interval.

### primary trait summary statistics file
Primary trait file is named by *TraitName*.*CHR*.*STARTBP*.*ENDBP*.txt. It is space-delimited and has CHR, BP, and SNP. It also has to carry one of STAT, T, or P. If P is specified, the two-sided P-value will be transformed into Z-statistic. STAT or T will be approximated as Z-statistic.

### reference LD file
One reference LD files should be provided for each interval. It is a tab-delimited file without a header. The file name is specified as locus.*CHR*.*STARTBP*.*ENDBP*.txt.gz. Each row is a marker, and it contains the following columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, and is followed by two alleles for each individual. 

## Run JLIM 2.0 permutation module on the example data
`CommandsExample.sh` will execute the permutation module sequentially.

```bash 
cd examples_2.0/

./CommandsExample.sh
```

This will create the following intermediate files for JLIM: 
- data file
- snps file
- positions file
- dosage file
- assoc file
- vars file
- betas file
- permutation file
- meta.assoc file

## Running JLIM on your data

First, assemble all loci you want to analyze and create one peaks file per chromosome, with a comma-separated list of all the genomic coordinates of the focal SNPs. For each cohort obtain a samples file and for each cohort/chromosome, obtain a corresponding bimbam and map file. Make sure you prepare all files mentioned the above. 

### 1) Generate the reference LD
*At this point, we recommend users to migrate to downloadable reference genotype panels provided with JLIM 2.5. However, JLIM 2.5 still supports the reference LD files in previously versions prepared as follows.* We provide a sample script to extract LD info for non-Finnish Europeans from downloaded 1000 genomes project vcf files. For example, if the vcf files are present in /data/1000genomes/ftp/release/20130502.

```bash
fetch.refld0.EUR.pl /data/1000genomes/ftp/release/20130502/ primary_trait/indexSNP.tsv ld0/
```

### 2) For each cohort/chromosome, run the “CutBimbam.py” script:

```
python CutBimbam.py bimbam_file map_file peaks_file samples_file output_string bimbam_separator_string map_separator_string map_position_int Maf_float interval_int
```

With the corresponding arguments:
- **output_string** cohort_name.chr, where chr is the chromosome number. This will be the beginning of the name of the output files: the snps and data files. Use the name of the cohort and the chromosome number separated by a “.” as in the example. You can include the folder where these files will be stored as a prefix (see example).
- **bimbam_separator_string** This is the character that separates columns in the bimbam file.
- **map_separator_string** This is the character that separates columns in the map (or info) file.
- **map_position_int** In the map file, the column number that holds the genomic positions will be denoted by map_position_int. The first column corresponds to a 0.

Optional arguments:
- **Maf_float** Default 0.05. This is the minor allele frequency (MAF) cutoff . The script will eliminate all SNPS with a minor allele frequency below this cutoff. 
- **interval_int** Default 200,000. This is the size of the interval (in BP), which will be analyzed by JLIM

Output:  
The script will produce a ‘data’ and a ‘snps’ file for each peak in the chromosome. These will be named with the output_string and the chromosome coordinates as: **output_string.start_BP.end_BP.snps.gz**, where start_BP and end_BP are the chromosome coordinates of the boundaries of the interval analyzed. These will in turn be the coordinates of each focal SNP listed in the peaks file +/- interval_int/2.

### 3) Make directories
For each locus, make  a directory where all files with meta-analyzed statistics will be stored, so that JLIM can use the directory to obtain the permutation p-values. These folders should follow the naming format locus.chr.startBP.endBP. We recommend storing these folders in a directory with the secondary phenotype’s name.

```bash
mkdir secondary_phenotype_name
mkdir secondary_phenotype_name/locus.chr.startBP.endBP
```

### 4) For each locus, run the “Makedosage.py” script:
This script will merge all genomic data in the different ‘data’ and ‘snps’ files to make a dosage file with the genotypes of all samples in the locus. It also makes locus-cohort specific position files.
```bash 
python Makedosage.py cohort_number_int chr  secondary_phenotype_name/locus.chr.startBP.endBP/dosage_file cohortA_data_file cohortA_snps_file cohortB_data_file cohortB_snps_file
```

With the corresponding arguments:
- **cohort_number_int** The number of cohorts that will be meta-analyzed together
- **chr** The chromosome number
- **dosage_file** This file will have the genotypes of all samples in the locus, and will be used to calculate in-sample LD by JLIM. Make sure to include the appropriate path to where the file should be stored
- **cohortA_data_file** and **cohortA_snps_file** list all data and snps files from the cohorts that will be meta analyzed together. The script expects cohort_number_int data files and cohort_number_int snps files

### 5) Run regressions
For each locus-cohort-secondary phenotype, run the “RunRegressions.py” script:  
This script will obtain summary statistics for each locus. It will generate permutation files with estimated effect sizes (betas_files) and estimated variances (vars_files), as well as individual summary statistics files (assoc_files). 

``` bash
python RunRegressions.py data_file snps_file samples_file secondary_phenotype_file covariates_file regression_id_string covariates_list permutation_number_int chr cohort_number_int
```

With the corresponding arguments:
- **data_file** and **snps_file** These files should be locus and cohort specific 
- **samples_file**, **phenotypes_file**, and **covariates_file** These files should be cohort specific
- **regression_id_string** This argument will be used in the naming of the outputs of this script. It should include the cohort name, the locus ID (chr.start.end) and the name of the secondary phenotype included in the secondary phenotype file. It can also include a folder where the outputs will be stored. Ideally: cohort_name.chr.start_bp.end_bp.secondary_phenotype_name
- **covariates_list** This is a comma-separated list, which indicates the covariates from the covariate_file that will be used in the regression. The left-most column (column 0) should correspond to the sample names. The first covariate column after the sample name would be column 1.
- **permutation_number_int** The number of permutations that will be generated. This determines the lower bound and precision of the JLIM p-value.
- **chr** The chromosome number
- **cohort_number_int** The number of cohorts that will be meta-analyzed together

### 6) Meta analyze
For each locus-secondary phenotype, when all of the cohort specific regressions have been finished, run “METAmergecohorts.py” to combine the cohort specific statistics into meta-analyzed files. These will include a summary statistic association file (meta.assoc_file) and a permutation file (meta.dump.all)

```bash
python METAmergecohorts.py cohort_number_int assoc_file secondary_phenotype_name/locus.chr.startBP.endBP/meta_id_string beta_file_1 vars_file_1 betas_file_2 vars_file_2 ...

```

For example if two cohorts (A and B) are being meta analyzed:

```
python METAmergecohorts.py 2 CohortA_assoc_file secondary_phenotype_name/locus.chr.startBP.endBP/meta_id_string CohortA_beta_file CohortA_vars_file CohortB_beta_file CohortB_vars_file
```

with the following arguments:
- **cohort_number_int** The number of cohorts that will be meta-analyzed together
- **CohortA_assoc_file** Use any assoc_file from the locus. The specific cohort that generated this file should be irrelevant.
- **meta_id_string** This argument will be used in the naming of the outputs of this script. It should include the set of cohorts combined, the locus ID (chr.startBP.endBP) and the name of the secondary phenotype included in the phenotype file. Make sure to include the appropriate path to where the files should be stored
- **CohortA_beta_file** and **CohortA_vars_file** The beta and vars files generated from the Runregression.py script. Each beta file should proceed the var file from the same cohort. Pairs of beta and vars files from all cohorts that will be meta-analyzed should be provided.

### 7) Run JLIM
Finally, for each config file created (one per primary-secondary phenotype pair), run JLIM v2.5 as below:

```bash
run_jlim.sh --index-snp 1:1071318  --maintr-file Height/Height.1.1000000.1200000.txt \
--sectr-file Phen1/locus.1.1000000.1200000/CohsABC.1.1000000.1200000.Phen1.meta.assoc.linear.gz \
--ref-ld ld0/locus.1.1000000.1200000.txt.gz --sectr-ld Phen1/locus.1.1000000.1200000/CohsABC.dosage.gz \
--perm-file Phen1/locus.1.1000000.1200000/CohsABC.1.1000000.1200000.Phen1.meta.mperm.dump.all.gz \
--manual-window-boundary 1000000-1200000 --output-file 1.1000000.1200000.out

run_jlim.sh --index-snp 1:1268710 --maintr-file Height/Height.1.1200000.1400000.txt \
--sectr-file Phen1/locus.1.1200000.1400000/CohsABC.1.1200000.1400000.Phen1.meta.assoc.linear.gz \
--ref-ld ld0/locus.1.1200000.1400000.txt.gz --sectr-ld Phen1/locus.1.1200000.1400000/CohsABC.dosage.gz \
--perm-file Phen1/locus.1.1200000.1400000/CohsABC.1.1200000.1400000.Phen1.meta.mperm.dump.all.gz \
--manual-window-boundary 1200000-1400000 --output-file 1.1200000.1400000.out
```

### Excluding samples
If you need to exclude any samples, you should exclude them from the secondary phenotype file, or set their phenotype to “nan”. This file can only take quantitative values (or nan). The sample file should match the genotypes in the bimbam file exactly, so no samples can be excluded from this file.
