# selfie
all script associated with gf-selfie

### VCF cleaning

This script cleans a given VCF file by removal of low confidence variants. The following criteria are applied:

1) Variants are kept if they carry the 'PASS' flag assigned by GATK.

2) The maximal length of InDels and MNP is set to 1kb.

3) The minimal number of reads at a variant position is 20 and the maximum number of reads is 500. This is a coverage cutoff.

4) The variant is biallelic.


```
python vcf_cleaner.py
--in <INPUT_VCF>
--out <OUTPUT_VCF>
```
`--in` species the VCF input file.

`--out` species the VCF output file.


### Selection of heterozygous variants

This script extract heterozygous variants from a given VCF file. Variants are considered heterozygous if the alternative allele frequency is between 0.1 and 0.9. 

```
python select_hetero_SNVs.py
--in <VCF(INPUT)>
--out <VCF(OUTPUT)>
```
`--in` species the VCF input file.

`--out` species the VCF output file.


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



### Plotting significant variant frequencies

This scripts generates a plot for a pseudochromosome or a specified target region. The number of all variants, the number of significant variants and the proportion of significant variants are displayed.

```
python sig_var_around_candidate_region.py
--sig <VCF_WITH_SIGNIFICANT_VARIANTS>
--all <VCF_FILE_WITH_ALL_VARIANTS>
--seq <NAME_OF_SEQ_OF_INTEREST>
--fig <FIGURE_FILENAME>

optional:
--window <SIZE_OF_SLIDING_WINDOW>
--step <STEP_SIZE>
--x1 <INTERVALL_STAR_POSITION>
--x2 <INTERVALL_END_POSITION>
--score <ACTIVATES_P_VALUE_BASED_SCORE>
``` 

`--sig` species a VCF input file which contains all variants with a significant difference in the distribution of reads in both pools.

`--all` species a VCF input file with all variants.

`--seq` species the name of a pseudochromosome or sequence. Only variants on this sequence will be plotted.

`--fig` species the figure output file name. The file extension determines the file type. Frequently used options are PDF, PNG, JPG, and SVG.

`--window` species the size of a sliding window which is used to summarize variant numbers. Each window should contain many variants to achieve a good representation. Setting the window size to extremely small values can result in artifacts.

`--step` species the size of steps by which the sliding window is moved along a pseudochromosome/sequence.

`--x1` species the start position of a region of interest.

`--x2` species the end position of a region of interest.

`--score` this flag activate the inclusion of a adj. p-value-based score. The log10() of the average adjusted p-value of all variants in a window is displayed. 




### Script for extraction of high impact variants from SnpEff results

This script extracts a customized selection of high impact variants from the SnpEff results file. Reference: [10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321).

```
python SnpEff_result_parser.py
--in <FULL_PATH_TO_INPUT_VCF>
--out <FULL_PATH_TO_OUTPUT_TEXT_FILE>
--gff <GFF_FILE_FOR_GENE_ID_MAPPING>

optional:
--anno <ANNOTATION_FILE>
``` 


`--in` species a VCF file which was annotated by SnpEff.

`--out` species an output text file where the high impact variants will be summarized.

`--gff` species a GFF file for the mapping of gene IDs.

`--anno` species a text file containing functional annotations of the genes. The first column is expected to contain the gene IDs with functional annotations in the following columns.




### Script for selection of high impact variants in target region

This script selects high impact variants predicted by SnpEff in a specific region. The input file should only contain variants of one pseudochromosome (sequence).

```
python select_candidates.py
--anno <FUNCTIONAL_ANNOTATION_FILE>
--gff <GFF3_FILE>
--in <HIGH_IMPACT_SNPEFF_FILE>
--out <OUTPUT_FILE>
--start <START_POSITION>
--end <END_POSITION>
``` 

`--anno` species a text file with functional annotation of the genes in the GFF file.

`--gff` species a GFF file which contains information about the positions and structures of annotated genes.

`--in` species a text file which contains information about high impact variants.

`--out` species a output file which contains information about high impact variants.

`--start` species the start position of the interval of interest.

`--end` species the end position of the interval of interest.





### References
