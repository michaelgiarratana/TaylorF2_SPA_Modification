#!/bin/bash

python PE_From_Frames.py \
  --mass1 36.0 \
  --mass2 29.0 \
  --s1 0.5 \
  --s2 0.5 \
  --distance 400.0 \
  --inclination 0 \
  --eventT 1197008880 \
  --H_path /content/X-BBH1_END_H1-1197008864-16.gwf \
  --L_path /content/X-BBH1_END_L1-1197008864-16.gwf \
  --H_channel BBH1_END_H1 \
  --L_channel BBH1_END_L1 \
  --out_dir visualising_the_results \
  --run_label BBH1_END_M_eta_dL_a1_a2 \
  --source_type BBH
