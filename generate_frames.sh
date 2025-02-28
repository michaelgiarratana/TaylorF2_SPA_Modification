#!/bin/bash

python3 generate_signal.py \
	--config parameters.txt \
	--psd_file lalsim.SimNoisePSDaLIGOZeroDetHighPowerGWINC \
	--outputfile outputfile.gwf \
	--channel TestMOD \
	--det H1 \
	--eventT 1197008880 \
	--ra 0 \
	--dec 0 \
	--psi 0 \
	--srate 16384 \
	--pad 8
