#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains figure 2 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 2A - Coding mutation burden
Rscript $SRCDIR/plt-fig2a-tmb-summary-barplot.r \
     --in_file=$TMB_SUMMARY \
     --id_map=$ID_MAPPING \
     --sample_order=$FIG_2_ORDER \
     --out_file=$SRCDIR/fig2a-tmb-summary-barplot.pdf


## Figure 2B - Oncoprint 
## Edit plt-fig2b-oncoprint.html to include metadata.txt and oncomatrix.txt inputs

## Figure 2C - COSMIC SBS/DBS/ID contributions
## Edit plt-fig2c-mutational-signature-barplot.html file to include SBS, DBS, and ID signatures