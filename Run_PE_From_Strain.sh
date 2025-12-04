#!/bin/bash

N_RUNS=20

BASE_LABEL="BBH23_MULTISEED_phic_mc_q_a1_a2_MOD_0_CE"

# Fixed Seed: 88170235

SEEDS=(58541173 65165548 88595755 61509788 \
	58133785 39670717 53084943 83849832 \
	29058102 87926741 79000453 45285690 \
	88576689 76912072 30319968 66971490 \
	76559093 24418224 84875658 42796138)

for i in $(seq 1 $N_RUNS); do
    echo "======================================"
    echo " Starting Run $i / $N_RUNS "
    echo "======================================"

    RUN_LABEL="${BASE_LABEL}_${i}"
	SEED=${SEEDS[$((i-1))]}
	python PE_From_Strain.py \
		--mass1 35 \
		--mass2 15 \
		--s1 0.3 \
		--s2 0.3 \
		--distance 700.0 \
		--inclination 0 \
		--eventT 100 \
		--phi_ref 1.3 \
		--f_min 10.0 \
		--f_max 112.202 \
		--f_ref 0.0 \
		--out_dir /home/giarratana/TaylorF2_SPA_Modification/october_results/ \
		--run_label "${RUN_LABEL}" \
		--source_type BBH \
		--injection_type pickle \
		--detector CE \
		--strain_path /home/giarratana/TaylorF2_SPA_Modification/saved_dict_BBH23_MOD_0_CE.pkl \
		--duration 22 \
		--begin 11 \
		--seed ${SEED}
	echo "Run ${i} completed. Results saved/"
    echo
done

echo "All ${N_RUNS} runs completed successfully."
