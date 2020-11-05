# selfie
all script associated with gf-selfie

### Filtering of VCF file based on a gold standard

This scripts reduces the variants in a VCF file to those positions which are supported by a gold standard.

```
python filter_VCF_by_gold_standard.py
--vcf <VCF (INPUT)>
--gold <GOLD_VCF (INPUT)>
--out <VCF (OUTPUT)>
``` 

`--vcf` species the VCF input file. All variants in this file are checked against the gold standard VCF file. Only variants at matching positions are written into the output file.

`--gold` species the gold standard VCF file. All variants in the input file are checked against this. Only variants with support from this gold standard are written into the output file.

`--out` species the VCF output file. All variants of the input file which are supported by the gold standard are written into the output file.




### Delta allele frequency plots

This script is used to visualize the distribution of delta allele frequencies across pseudochromosomes.

```
python delta_allele_frequency.py\n
--input_vcf <FILENAME>
--reference_file <FILENAME>
--output_dir <DIRECTORY_NAME>[will be generated if required]
--pool1 <sample name in VCF; multiple samples names can be provided comma-seperated>
--pool2 <sample name in VCF; multiple samples names can be provided comma-seperated>
					
optional:
--minP1cov <FLOAT, lower coverage cutoff for pool1>
--maxP1cov <FLOAT, upper coverage cutoff for pool1>
--minP2cov <FLOAT, lower coverage cutoff for pool2>
--maxP2cov <FLOAT, upper coverage cutoff for pool2>
``` 

`--input_vcf` species VCF file with two data columns, one for each of the pools.


`--reference_file` species a reference FASTA file.


`--output_dir` species a output directory where all figure files will be placed. This folder will be created if it does not exist already.

`--pool1` species the name of pool1 in the VCF file. This string needs to match the header of the respective column in the VCF.

`--pool2` species the name of pool2 in the VCF file. This string needs to match the header of the respective column in the VCF.

`--minP1cov` species the minimal coverage (number of reads) of pool1.

`--maxP1cov` species the maximal coverage (number of reads) of pool1.

`--minP2cov` species the minimal coverage (number of reads) of pool2.

`--maxP2cov` species the maximal coverage (number of reads) of pool2.









### Fisher's exact test

fisher_exact_test_corrects_for_multiples_testing.py


### References
