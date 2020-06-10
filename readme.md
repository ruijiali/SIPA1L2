SIPA1L2 as a risk factor implicated in Alzheimer’s disease

Our investigation demonstrated the presence of loci of susceptibility for AD in SIPA1L2.

## 1.Pre-requisite software

Table 1. List of pre-requisite software and the available information.

| No.  | Software  | Version | Availability                                                 |
| :--- | --------- | ------- | ------------------------------------------------------------ |
| 1    | Michigan  | 1.2.4   | https://imputationserver.sph.umich.edu/index.html            |
| 2    | PLINK     | 1.9     | http://www.cog-genomics.org/plink2/                          |
| 3    | LiftOver  |         | https://genome.ucsc.edu/cgi-bin/hgLiftOver                   |
| 4    | SHAPEIT   | 2.r900  | https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html |
| 5    | HRC check | 4.2.11  | https://www.well.ox.ac.uk/~wrayner/tools/                    |
| 6    | BCFTOOLS  | 1.9     | http://samtools.github.io/bcftools/bcftools.html             |
| 7    | ANNOVAR   |         | http://annovar.openbioinformatics.org/en/latest/             |

## 2.Methods

### 2.1 pre-processing steps

1. Pre-QC executed by PLINK.

   ```
   $Plink --bfile file1 --keep-allele-order --hwe 0.00001 --geno 0.05 --maf 0.005 --max-maf 0.01 --mind 0.05 --make-bed --out file2
   ```

2. Use LiftOver to convert the genome coordinates into hg19.

   ```python
   python liftOverPlink.py -m file.map -p file.ped -o plink_ped -c hg18ToHg19.over.chain.gz -e liftOver
   ```

   


### 2.2 Genotype imputation

Michigan Imputation service (https://imputationserver.sph.umich.edu/index.html) .

### 2.3 Post-processing steps and QC.

1. Filter SNPs with R2 < 0.8.

   ```
   for i in {1..22}; do bcftools view -i 'R2>.8' -Oz chr${i}.dose.vcf.gz > file1; done
   ```

2. Changing the ID to '%CHROM:%POS:%REF:%ALT'.

   ```
   for i in {1..22}; do bcftools norm -Ou -m+any file1| bcftools annotate --output-type b --output file2  -I '%CHROM:%POS:%REF:%ALT'; done
   ```

3. Convert bcf format to plink format.

   ```
   $PLINK --bcf file2--keep-allele-order --allow-extra-chr 0 --const-fid --split-x b37 no-fail --vcf-idspace-to _ --make-bed --out file3 
   ```

6. Merge chromosomes.

   ```
   $PLINK --bfile file3 --merge-list mergelist.txt --biallelic-only --make-bed --out file4
   ```

5. Post-QC executed by PLINK.

   ```
   $PLINK --bfile file4 --keep-allele-order --hwe 0.00001 --geno 0.05 --maf 0.005 --max-maf 0.01 --make-bed --mind 0.05 --noweb --out file5
   ```

   


### 2.4 Gene-based functional annotation.

1. Convert plink to VCF4 file.

   ```
   $PLINK --bfile file6 --recode vcf-iid --out file7
   ```

2. ANNOVAR annotation

   ```
   perl convert2annovar.pl -format vcf4old ./Michi_filter_final_dup.vcf > ./Michi_filter_final_dup.avinput
   
   perl annotate_variation.pl -geneanno -dbtype refGene -buildver hg19 ./Michi_filter_final_dup.avinput ./humandb/
   ```

### 2.5 Statistical analysis

"SIPA1L2.R" can be downloaded  from https://github.com/ruijiali/SIPA1L2. 

### 2.6 Bioinformatic analyses

Table 2. List of bioinformatics analysis tools and the available information.

| No.  | Tools                                      | Availability                                          |
| :--- | ------------------------------------------ | ----------------------------------------------------- |
| 1    | NCBI                                       | https://www.ncbi.nlm.nih.gov/gene/                    |
| 2    | GeneCards database                         | https://www.genecards.org/                            |
| 3    | targetvalidation                           | https://www.targetvalidation.org/                     |
| 4    | varsome                                    | https://varsome.com/                                  |
| 5    | genepine                                   | http://grch37.genepipe.ncgm.sinica.edu.tw/variowatch/ |
| 6    | Noncode                                    | http://www.noncode.org/                               |
| 7    | hosphosite                                 | https://www.phosphosite.org/homeAction                |
| 8    | gnomad                                     | https://gnomad.broadinstitute.org/                    |
| 9    | UCSC Genome Browser on Human               | https://genome.ucsc.edu/                              |
| 10   | Online Mendelian Inheritance in man (OMIM) | https://omim.org/entry/                               |