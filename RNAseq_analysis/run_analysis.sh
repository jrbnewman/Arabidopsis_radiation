#!/bin/sh
sas -work /mnt/store/sas_work -memsize 20g -sysin ./analysis_arabidopsis_fusions_run_model_trt_time.sas
sas -work /mnt/store/sas_work -memsize 20g -sysin ./analysis_arabidopsis_fusions_run_model_trt_time_main_effects.sas
sas -work /mnt/store/sas_work -memsize 20g -sysin ./analysis_arabidopsis_fusions_run_model_trt_time_contrasts.sas
sas -work /mnt/store/sas_work -memsize 20g -sysin ./analysis_arabidopsis_up_down_sign_by_contrast_fusion.ssas
sas -work /mnt/store/sas_work -memsize 20g -sysin ./export_nominal_summaries.sas



