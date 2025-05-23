#!/bin/bash

nohup python3 main_PE.py --psd-dir=/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/3dkde_studies/psd/  \
				   --mass1=1.588621354587203   --mass2=1.378895267390714   --luminosityDistance=282   \
				   --eos=APR4_EPP  --waveform-approximant=IMRPhenomPv2_NRTidal  --search-waveform-approximant=TaylorF2  \
				   --ra=0.9690181423973839 --dec=1.0485191476021183 --inc=2.56171991199752 --psi=0.3643211321014299 \
				   --chi1=-0.008104654230022925 --chi2=0.004063180021609744 > out1.dat& 

nohup python3 main_PE.py --psd-dir=/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/3dkde_studies/psd/  \
				   --mass1=1.588621354587203   --mass2=1.378895267390714   --luminosityDistance=100   \
				   --eos=APR4_EPP  --waveform-approximant=IMRPhenomPv2_NRTidal  --search-waveform-approximant=TaylorF2  \
				   --ra=0.9690181423973839 --dec=1.0485191476021183 --inc=2.56171991199752 --psi=0.3643211321014299 \
				   --chi1=-0.008104654230022925 --chi2=0.004063180021609744 > out2.dat& 
