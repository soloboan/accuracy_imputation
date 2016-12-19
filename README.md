<img src="icon.png" align="right" />
## Computing accuracy of imputation
Computing imputation accuracy with imputed data formats ranging from 

- dosages (discrete genotypes 0, 1, 2)
- gene content (genotypes ranging from 0-2)
- allele format (mostly in linkage and transposed linkage formats eg. PLINK PED and TPED)

Imputation accuracy is computed using four methods 
- correlation between imputed and true genotypes
- Allelic correct rate (1 - allelic error rate)
- Percentage of correct calls
- correlation between centered-imputed and centered-true genotypes (see Calus et al. 2014 - Evaluation of measures of correctness of genotype imputation in the context of genomic prediction: a review of livestock applications). 

# The following are the R-script argument 
Truegeno == filename of true genotype (PLINK or allele dosage or gene content format).

Imputedgeno == filename of true genotype (PLINK or allele dosage or gene contetn format).
  	-note that both imputedgeno and truegeno should have the same file format.

mapHD == filename of the HIGH density SNP map file.

mapLD == filename of the LOW density or imputed genotype SNP map file.

nIID == number of animals genotyped, you can specify higher number but not less. 

format == The format of the genotype data (either PLINK or gene content/allele dosage format).

accuracytype == The method to use in calculating accuracy (1. simple 2. calus_mulder 3. both). Read: Evaluation of measures of correctness of genotype imputation in the context of genomic prediction: a review of livestock applications (Calus et al. 2014 - http://dx.doi.org/10.1017/S1751731114001803). 

missgeno == the value of the missing genotype in the data file (eg. NA, 5, -9 ...).


