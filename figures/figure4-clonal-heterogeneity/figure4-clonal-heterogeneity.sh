#!/bin/bash
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper

################################################################# /COPYRIGHT ###
################################################################################

## This bash script contains figure 4 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 4A - MutationTimer summary barplots
Rscript $SRCDIR/plt-fig4a-mutation-timer-summary-barplot.r \
  --in_file=$MTSUMMARY \
  --pp=$PURITY_PLOIDY \
  --id_map=$ID_MAPPING \
  --out_file=$SRCDIR/ig4a-mutation-timer-summary-barplot.pdf


## Figure 4B - WCM12 alluvial plot
Rscript $SRCDIR/plt-fig4b-wcm12-alluvial-plot.r \
  --vcf=$MTDIR/PM12-Z10-1-Case-WGS--PM12-EBC2-1-Ctrl-WGS-mutationtimer-snv.vcf,$MTDIR/PM12_Z4_2_Case--PM12_EBC2_2_Ctrl-mutationtimer-snv.vcf,$MTDIR/PM12_Z13_1_Case--PM12_EBC2_2_Ctrl-mutationtimer-snv.vcf \
  --vcf_name=WCM12_PR_1,WCM12_BR_3,WCM12_BR_2 \
  --out_file_pdf=$SRCDIR/fig4b-wcm12-alluvial-plot.pdf
