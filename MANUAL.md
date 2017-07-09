## RFMIX v2.XX - Local Ancestry and Admixture analysis

### About

RFMIX is a program to identify the ancestry of genomic segments using random forest discriminative machine learning methods combined with a conditional random field model of the linear chromosome. This method is capable of determining both the dominant single ancestry of non-admixed individuals as well as decomposing the genome of admixed individuals into discrete segments of arbitrary size originating from different ancestral populations. The ancestral populations are represented to the program as a reference panel of present day populations thought to be direct decedents of the populations which contributed to today's admixed populations.

### LICENSE

RFMIX v2.XX is available for use without fee for academic research use only. For users wishing to use RFMIX for commercial purposes, explicit licenses may be granted for such use. Please contact [Stanford Office of Technology Licensing](http://otl.stanford.edu/) for more information.

### Requirements

+ A **phased** VCF/BCF file containing "query" haplotypes which are to be analyzed.
+ A **phased** VCF/BCF file containing reference haplotypes (in any order)
+ A reference sample map file matching reference samples to their respective reference populations
+ A genetic map file
+ All files mapped or referenced to the same genome assembly
+ All diploid genotypes must be phase resolved
+ In addition to the above, you must at minimum also specify a basename (prefix) for output files, and the chromosome to analyze from the VCF/BCF inputs, even if they contain only one chromosome
+ If BCF files are used, *bcftools* must be installed and available in the PATH environment setting
+ VCF/BCF files may be gzip compressed, and should be indexed using *bcftools*

### QUICKSTART

Run the program with no command line options to see a list of options accepted and a terse description of what they do. This help message is intended primarily as a reminder for how to run the program.

The following options are required:

~~~~~~~~~~~~
	-f <query VCF/BCF file>
	-r <reference VCF/BCF file>
	-m <sample map file>
	-g <genetic map file>
	-o <output basename>
	--chromosome=<chromosome to analyze>
~~~~~~~~~~~~

RFMIX will automatically perform an intersection of the query VCF and the reference VCF to determine SNPs which are in common between them. You do not need to do this yourself or do this for each specific combination of a query VCF file and various reference panels you might want to use. The 6 required options above are all that is needed to run the program.

It is recommended that BCF files be used as input. The samples in the VCF/BCF files may appear in any order. It is recommended that the entire genome (all chromosomes) be contained in one VCF/BCF file for the query, and the reference rather than seperate by chromosome. If the BCF files are indexed, RFMIX will skip directly to the chromosome being analyzed.

The genetic map file is tab delimited text containing at least 3 columns. The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM. Any number of columns or other information may follow, it is ignored. The chromosome column is a string token (which may be an string of digits) that must match those used in the VCF/BCF inputs. The genetic map file should contain the map for the entire genome (all chromosomes). Blank lines and lines beginning with a '#' are ignored.

The sample map file specifies which subpopulation each reference sample represents. It is tab delimited text with at least two columns. The first column gives the sample name or identifier, which must match the one used in the reference VCF/BCF. The second column is a string naming a subpopulation and may contain spaces (e.g., "European", or "East African"). RFMIX will assign all distinct subpopulation names it finds in the sample map file an index number, in alphabetical order. The output will reference by index number; the order is given at the top of the output files. Blank lines and lines beginning with a '#' are ignored in the sample map file. Prefixing a sample with either # or \^ will exclude the sample from the reference input without needing to remove it from the reference VCF/BCF. Any sample not defined in the sample map will not be loaded from the reference VCF/BCF. This is a simple way to manipulate the content of your reference data and include or exclude entire subpopulations.

### Output files

RFMIX upon completion will output two main files of interest: the most likely assignment of subpopulations per CRF point (\<output basename\>.msp.tsv), and the marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point (\<output basename\>.fb.tsv). The latter is produced by computing the forward-backward algorithm on the CRF, and the former by using the Viterbi algorithm. The .msp.tsv file is condensed such that CRF windows are combined if all query samples are in the sample subpopulations for successive windows. Thus, each line might represent several CRF points.

Both output files are tab separated values forming a matrix with rows corresponding to genomic position and columns corresponding to haplotypes. The files include column headers and leading columns that indicate the position or range covered for each row. For the forward-backward results, haplotypes are tab delimited, but the array of probabilities for each haplotype at each window (row) is a set of space delimited columns within each tab delimited haplotype column. The order and names of the reference subpopulations are indicated in a header row.

Global *diploid* ancestry estimates are computed by RFMIX and output to \<output basename\>.rfmix.Q, corresponding to the .Q output files of the global ancestry analysis programs *fastStructure* or *ADMIXTURE*. Please note that the ordering of subpopulations in columns for *fastStructure* or *ADMIXTURE* is not guaranteed to be the same as RFMIX and the ordering may change for different runs of those programs. RFMIX will always output in the same order (indicated in a header line) if the reference file used is the same. Additionally, RFMIX will add the sample name/id from the VCF input as a leading column.

If EM is used (see below), output is generated on completion of each iteration, but each successive iteration overwrites the results of the previous. At completion of the program, the output files contain the final results from the last iteration.

### Further options of interest

Additional options of interest are the CRF spacing size, or the number of SNPs each point of conditional random field model represents (-c \<# of SNPs\>), and the random forest window size (-r \<# of SNPs\>). Either of these options may be specified instead as a genetic distance in cM. If the value is less than 1.0 (2.0 for -r), it is interpreted as a genetic distance. Otherwise, it is interpreted as the number of SNPs. The CRF spacing size must be less than or equal to the random forest window size, the program will automatically expand random forest window sizes to include any SNPs in the input that would fall between windows otherwise. These parameters have default values and do not need to be specified, but it is generally desired to control this explicitly.

The number of generations since contact of populations composing the query or putatively admixed samples can be specified by -G \<# of generations\>. This parameter is not required and will default to 8 generations.

**IMPORTANT:** If the CRF window size and random forest window sizes are not equal, an additional parameter is needed for accurate results. At the time of writing, this parameter is not automatically learned by the program and thus must be specified manually on the command line (-w \<decimal value\>). The value must be greater than zero. Future versions of this program will perform an internal simulation to automatically determine an optimal value. At present, this is not done. If not specified explicitly, the program will set the value to 2.0 * RF size / CRF size. The optimal value may be dependent on the content and SNP density of the reference panel.

RFMIX can iterate the algorithm to perform expectation-maximization (EM) optimization of the model, by including the query/admixed and their current reference subpopulation probabilities, per each CRF window, as additional reference haplotypes starting with the first EM iteration after the initial analysis. This allows novel haplotypes found in the query data to inform classification and improve results. This may make substantial improvements if the reference populations used are proxies for ancestral populations which are not available to be used as reference directly. This option is off by default. To turn EM on, use -e \<max # of iterations\>. EM will stop when the specified number of iterations is reached, or the results no longer change from one iteration to the next, whichever happens first.

In the case a set of reference haplotypes may not be of "pure" ancestry and may themselves be somewhat admixed, the option --reanalyze-reference will cause the program to analyze the reference haplotypes as if they were query haplotypes, in addition to analyzing the query input. On the first EM iteration and further, the results of this analysis will update and replace the uniform assignment of the entire chromosome to a single reference population. Thus, segments of some reference haplotypes may be reassigned to a different subpopulation than the one specified in the reference sample map (-m). To use this, both --reanalyze-reference and -e \<# of iterations\> options must be specified. --reanalyze-reference has no effect if EM is not turned on.

The option --analyze-range=\<string\> can be used to restrict analysis only to the range of positions given in Mbp. For instance --analyze-range=\<30.5-50\> will analyze only the portion of the chromosome falling within the range 35,000,000 to 50,000,000 bp. This may be useful to speed up the total analysis time when exploring how to use the program and its various options.

### Limitations

+ The quality and accuracy of the results depends directly on the extent to which the haplotypes provided for each reference population captures the breadth of genetic diversity within that population.
+ Single or small sample sizes often lack reference homozygote genotype calls at known variable sites and it is not possible to know from a VCF/BCF file whether read data supporting a reference homozygote call was observed. Unfortunately, this information is very informative to RFMIX and so it being censored may strongly impact results.
+ Accuracy of the results, and continuity of subpopulation assignment along each haplotype, depends on the accuracy of the phasing of diploid data which must be performed separately as a pre-requisite. The phase-correction features of RFMIX, if enabled, can not correct all phasing errors and will perform poorly if the initial phasing presented to the program is poor.
+ RFMIX version 2 accepts inputs containing missing data. However, at time of writing, this feature has not been tested or studied. The rate at which results degrade for increasing amounts of missing data is not known.
+ The program will assume genotypes are phased even if the VCF file indicates they are not.

### DISCLAIMER

For academic users or any other users using the program free of charge:

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM (RFMIX) "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
