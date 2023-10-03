libname wgbsloc '/TB14/TB14/sandbox/wgbs_cold_sandbox/sas_data';

data meth_data;
  set wgbsloc.meth_data2region;
  keep gene_id_cat region_type_cat strand_cat region_start_cat region_stop_cat chr pos;
run;

proc sort data=meth_data nodup;
by _all_;
run;

data wgbsloc.meth_pos2region;
  set meth_data;
  run;
  
