#!/bin/bash

python PE_From_Strain.py \
  --mass1 35 \
  --mass2 20 \
  --s1 0.5 \
  --s2 0.5 \
  --distance 400.0 \
  --inclination 0 \
  --eventT 1197008880 \
  --phi_ref 1.3 \
  --f_min 15.0 \
  --f_max 166.187 \
  --f_ref 0.0 \
  --out_dir /home/giarratana/TaylorF2_SPA_Modification/results/ \
  --run_label BBH4_tc_phic_Mc_q_a1_a2_MOD \
  --source_type BBH \
  --injection_type pickle \
  --strain_path /home/giarratana/TaylorF2_SPA_Modification/saved_dict_BBH4_MOD.pkl
