module load r

cd /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0006_WRAIR_kenread_ken

/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/MIPhackR_v01.R \
  --guide /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0006_WRAIR_kenread_ken/s0006_WRAIR_kenread_ken.xlsx \
  --target /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/targets_k13updated2023.tsv \
  --meta_site_key /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0006_WRAIR_kenread_ken/s0006_WRAIR_kenread_ken_site_key.tsv \
  --output /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0006_WRAIR_kenread_ken/output_extract_tsvs \
  --coverage_threshold 1 \
  --mutant_umi 1 \
  --substudy untreatedextracted \
  --sex both \
  --pre_treatment Y \
  --minimum_age ""\
  --maximum_age ""
  
  python3 /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/geoff_validation/GEOFF_tools_v02.py tsv_validate \
--yaml /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/geoff_validation/GEOFF_tools_parameters.yaml \
--study_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0006_WRAIR_kenread_ken/output_extract_tsvs/s0006_WRAIR_kenread_ken_1_study_data.tsv \
--validation_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/SAVED_validate_parameters_2024_09_24.tsv
