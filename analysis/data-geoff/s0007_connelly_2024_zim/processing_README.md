# 1. Load R

Before running the script, make sure to load the required version of R (4.0.0). You can load it using the following command:

```
ml r/4.4.0
```

# 2. Run the MIPhackR.R executable

To execute the MIPhackR.R script, use the following command, adjusting the paths and parameters as necessary for your study:
```
/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/MIPhackR_v01.R \
-g "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0007_connelly_2024_zim/input_mip/s0007_connelly_2024_zim.xlsx" \
-t "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/MIPhackR/targets_k13updated2023.tsv" \
-o "/nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0007_connelly_2024_zim/output_extract_tsv" \
-mut 1 \
-c 1 \
-m "" \
-s "untreatedextracted" \
-min_a "" \
-max_a "" \
-sex "both" \
-pre_t "Y"
```
Parameters:

-g (--guide, required): Path to the guide file generated from the study overview template. All of the meta fields at the end of the study_overview template should be filled out.

-t (--target, required): Path to the targeted mutation file. This file specifies the mutation targets. Do not modify this file (targets_k13updated2023.tsv).

-o (--output, required): Output directory where results will be saved. This should follow the structure: /GEOFF_META/studies_internal/sXXXX_lastname_year_studyqualifier.

-c (--coverage,required): The UMI count and coverage threshold. Set the appropriate numeric threshold for your analysis.

-m (--meta_site_key, optional): Metadata site key. Use this option if your metadata is stored in a file that is not attached to your sample_uids.

# 2. Run validation script

```
python3 /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/GEOFF_tools_v02.py tsv_validate \
--yaml /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/GEOFF_tools_parameters.yaml \
--study_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/studies/s0007_connelly_2024_zim/output_extract_tsv/s0007_connelly_2024_zim_v01_study_data.tsv \
--validation_tsv /nfs/jbailey5/baileyweb/bailey_share/GEOFF_META/GEOFF_code/SAVED_validate_parameters_2024_09_24.tsv
```