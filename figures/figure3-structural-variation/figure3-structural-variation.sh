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

## This bash script contains figure 3 associated scripts
set -euo pipefail
SRCDIR=$(realpath $(dirname $0))


## Figure 3A - FGA barplot
Rscript $SRCDIR/plt-fig3a-fga-barplot.r \
  --fga=$FGA \
  --tn_file=$TNFILE \
  --id_map=$ID_MAPPING \
  --site_colors=$SITE_COLORS \
  --purity_ploidy=$PURITY_PLOIDY \
  --out_file=$SRCDIR/fig3a-fga-barplot.pdf


## Figure 3B - Junction burden heatmap
Rscript $SRCDIR/plt-fig3b-junction-burden-heatmap.r \
       --in_file=$JUNCTION_SUMMARY_TABLE \
       --id_map=$ID_MAPPING \
       --patient_map=$PATIENT_COLORS \
       --site_map=$SITE_COLORS \
       --pp=$PURITY_PLOIDY \
       --tp53_status=$SOMATIC_ALTERATION_SUMMARY \
       --tmprss2_status=$TMPRSS2_STATUS \
       --sample_order=$SRCDIR/fig3a-fga-barplot.txt \
       --out_file=$SRCDIR/fig3b-junction-burden-heatmap.pdf


## Figure 3C - AR rearrangements (WCM90)
Rscript $SRCDIR/jabba/plotting/plt-fig3cd-jabba-genome-graphs.r \
    --samples=PM90-Z19-1-Case-WGS,PM90-Z21-1-Case,PM90-Z25-1-Case-WGS \
    --id_map=$ID_MAPPING \
    --jabba_gg=$JABBA_DIR/PM90-Z19-1-Case-WGS--PM90-EBC2-1-Ctrl-WGS.events.rds,$JABBA_DIR/PM90-Z21-1-Case--PM90-EBC2-1-Ctrl-WGS.events.rds,$JABBA_DIR/PM90-Z25-1-Case-WGS--PM90-EBC2-1-Ctrl-WGS.events.rds \
    --bed=$AR_BED \
    --padding=5E5 \
    --out_file=$SRCDIR/fig3c-AR-genome-graph.svg


## Figure 3D - RB1 rearrangements (WCM12)
Rscript $SRCDIR/jabba/plotting/plt-fig3cd-jabba-genome-graphs.r \
    --samples=PM12-Z10-1-Case-WGS,PM12_Z13_1_Case,PM12_Z4_2_Case \
    --id_map=$ID_MAPPING \
    --jabba_gg=$JABBA_DIR/PM12-Z10-1-Case-WGS--PM12-EBC2-1-Ctrl-WGS.events.rds,$JABBA_DIR/PM12_Z13_1_Case--PM12_EBC2_2_Ctrl.events.rds,$JABBA_DIR/PM12_Z4_2_Case--PM12_EBC2_2_Ctrl.events.rds \
    --bed=$RB1_BED \
    --padding=5E5 \
    --out_file=$SRCDIR/fig3d-RB1-genome-graph.svg

