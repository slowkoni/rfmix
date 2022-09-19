## RFMIX version 2

## LICENSE

NOTICE: **This software is available for use free of charge for academic research use only.** Commercial users, for profit companies or consultants, and non-profit institutions not qualifying as "academic research" must contact the [Stanford Office of Technology Licensing](http://otl.stanford.edu/) for a separate license to use RFMIX. This applies to this repository directly and any other repository that includes RFMIX source, executables, or git commands that pull/clone this repository as part of its function. Such repositories, whether mine or others, must include this notice. Academic users may fork this repository and modify and improve RFMIX to suit their research needs, but also inherit these terms and must include a licensing notice to that effect.

## Newer methods of Interest

High Resolution Ancestry Deconvolution for Next Generation Genomic
Data - [G-Nomix](https://github.com/AI-sandbox/gnomix)

## Changes made by Embark

The forward-backward result output is a tab separated file with the name \<output basename\>.fb.tsv

The first line is a comment line, that specifies the order of subpopulations:
eg:
```
#reference_panel_population: golden_retriever  labrador_retriever  poodle  poodle_small
```

The second line specifies the column names, and every following lines gives data on a chunk of the genome, called a conditional random field (CRF) point.

The first four columns specify the chromosome, genetic marker's physical position in basepair units and genetic position in centiMorgans, and the genetic marker's numerical index in the rfmix genetic map input file. The remaining columns give the probabilities that the CRF point for a genotype's haplotype was assigned to a specific reference panel population. A genotype has two haplotypes, so the number of probabilities for a genotype is 2*(number of reference panel populations). The number of columns in the file is 4 + (number of genotypes) * 2 * (number of reference panel populations.

## Building RFMIX

The quick way:
```sh
autoreconf --force --install # creates the configure script and all its dependencies
./configure                  # generates the Makefile
make
```

The long way (in case the quick way gives any trouble):
```sh
aclocal                      # creates aclocal.m4
autoheader                   # creates config.h.in
autoconf                     # creates configure
automake --add-missing       # creates Makefile.in
./configure                  # generates the Makefile
make
```

This will build `rfmix` and `simulate`.
