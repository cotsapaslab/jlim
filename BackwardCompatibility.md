# Backward compatibility to JLIM 1.0 and 2.0
Previous versions of JLIM allowed users to have full control of permutation process with a precomputed permutation file. We support this mode with the `--perm-file` option.

## Input file formats to run JLIM with precomputed permutation data

### Secondary trait LD file
This file should be a gzipped or plain-text genotype file in the Plink .ped format and with the identical set and order of SNPs as in the secondary trait association file. *Currently, JLIM does not allow missing genotypes in the secondary trait genotype file.* 

### Precomputed permutation file 
This is a gzipped Plink .mperm.dump.all file. This file should have the identical set and order of SNPs as in the secondary trait association file. The permutation files (e.g. `output_prefix.mperm.dump.all.gz`) can be generated as below: 

```bash
plink --bfile your_genotype_file --pheno your_phenotype_file --covar your_covarite_file \
--chr your_chr --from-bp your_start_bp --end-bp your_end_bp \ 
--linear --mperm 1000 --merpm-save-all --out output_prefix

gzip output_prefix.mperm.dump.all 
```

### Reference LD file
*At this point, we recommend users to migrate to downloadable reference genotype panels provided with JLIM 2.5.* However, JLIM 2.5 still supports the reference LD files prepared with `fetch.refld0.EUR.pl` in previously versions. This file is a tab-delimited file without column headings. Each row is a SNP. And it contains the 7 variant-information columns - CHROM, POS, ID, REF, ALT, QUAL, and FILTER, - and is followed by two columns for the alleles in each individual. For example, `fetch.refld0.EUR.pl` prepares the reference LD file for non-Finnish European population from the 1000 Genomes Project data as follows:

```
bin/fetch.refld0.EUR.pl /data/1000genomes/ftp/release/20130502/ examples/MS/MS-indexSNP.tsv examples/ld0/
```

for the genomic intervals defined in `examples/MS/MS-indexSNP.tsv` file, assuming that the 1000 Genomes Project data were downloaded into `/data/1000genomes/ftp/release/20130502/`. `fetch.refld0.EUR.pl` is a perl script and relies on **vcftools**. 

