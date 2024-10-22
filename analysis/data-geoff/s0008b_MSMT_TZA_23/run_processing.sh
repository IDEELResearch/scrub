#!/usr/bin/bash

ml r/4.4.0
/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/MIPhackR_v01.R \
-g "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0008b_MSMT_TZA_23/input_mip/s0008b_MSMT_TZA_2023.xlsx" \
-t "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/targets_k13updated2023.tsv" \
-o "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0008b_MSMT_TZA_23/output_extract_tsv" \
-mut 1 \
-c 1 \
-m "" \
-s "untreatedextracted" \
-min_a "" \
-max_a "" \
-sex "" \
-pre_t "Y"

python3 /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/GEOFF_tools_v02.py tsv_validate \
--yaml /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/GEOFF_tools_parameters.yaml \
--study_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0008b_MSMT_TZA_23/output_extract_tsv/s0008b_MSMT_TZA_23_v01_study_data.tsv \
--validation_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/SAVED_validate_parameters_2024_09_24.tsv
