import argparse
import sys
import math
import numpy as np
import lalsimulation as lalsim
import lal
from lalframe import FrWriteREAL8TimeSeries


# Argument parsing setup
parser = argparse.ArgumentParser(description="Process parameter configuration file.")
parser.add_argument("--config", required=True, help="Path to parameter configuration file.")
parser.add_argument("--psd_file", type=str, required=True, help="Path to the PSD file.")
parser.add_argument("--outputfile", required=True, help="Path to the output file.")
parser.add_argument("--channel", required=True, help="Channel name.")
parser.add_argument("--det", required=True, help="Detector name.")
parser.add_argument("--eventT", type=float, required=True, help="Event time.")
parser.add_argument("--ra", type=float, required=True, help="Right ascension.")
parser.add_argument("--dec", type=float, required=True, help="Declination.")
parser.add_argument("--psi", type=float, required=True, help="Polarization angle.")
parser.add_argument("--srate", type=float, required=True, help="Sample rate.")
parser.add_argument("--pad", type=float, required=True, help="Padding.")
args = parser.parse_args()  # Parse the arguments

# Initialize an empty dictionary for parameters
params = {}

with open(args.config, 'r') as f:  # Open the config file for reading
    for line in f:  # Loop through each line in the file
        if line.strip() and not line.startswith("#"):  # Skip lines starting w/ N 
            key, value = line.split(":")  # Split each line into key and value
            try:
                parsed_value = eval(value.strip())  # Convert value to its proper type
                
                # If eval returns a tuple with one element, extract the element
                if isinstance(parsed_value, tuple) and len(parsed_value) == 1:
                    parsed_value = parsed_value[0]
                    
                params[key.strip()] = parsed_value
            except Exception:
                params[key.strip()] = value.strip()  # Default to string if eval fails


channel=args.channel
det=args.det
eventT=float(args.eventT)
ra=float(args.ra)
dec=float(args.dec)
psi=float(args.psi)
srate=float(args.srate)
pad=float(args.pad)

hp, hc = lalsim.SimInspiralTDFromFD(**params)
hp.epoch += eventT
hc.epoch += eventT

detector = lalsim.DetectorPrefixToLALDetector(det)

h = lalsim.SimDetectorStrainREAL8TimeSeries(hp, hc, ra, dec, psi, detector) #pass detector object instead of string

t0 = int(h.epoch - pad)
seglen = math.ceil(h.deltaT * h.data.length + 2*pad) # Computes the event duration in seconds (w/ a buffer)

# Below is the procedure for injecting noise from a psd onto the hp, hc data
segdur = 16.0 # segment duration
srate *= 1
seglen = int(segdur * srate)
stride = seglen // 2  # stride between segments
start = segdur * np.floor(float(hp.epoch) / segdur)
end = segdur * np.ceil((float(hp.epoch) + hp.data.length * hp.deltaT) / segdur)
recdur = end - start # duration of data record
reclen = int(recdur * srate) # record length
numseg = 1 + (reclen - seglen) // stride # number of segments to make
data = lal.CreateREAL8TimeSeries(
    name = channel,
    epoch = start,
    f0 = 0.0,
    deltaT = 1.0 / srate,
    sampleUnits = lal.StrainUnit,
    length = reclen
)

seg = lal.CreateREAL8TimeSeries(
    name = 'noise',
    epoch = start,
    f0 = 0.0,
    deltaT = 1.0 / srate,
    sampleUnits = lal.StrainUnit,
    length = seglen
)

psd = lal.CreateREAL8FrequencySeries(
    name = 'psd',
    epoch = 0,
    f0 = 0.0,
    deltaF = 1.0 / segdur,
    sampleUnits = lal.SecondUnit,
    length = seglen // 2 + 1
)
  
# temp psd function to construct the psd calling command 
def psd_temp(func,i,ii):
  return func(i,ii)
psd_temp(eval(args.psd_file),psd,params['f_min'])


rng = lal.gsl_rng('mt19937', 123) # gsl random number generator

for i in range(numseg):
    lalsim.SimNoise(seg, 0 if i == 0 else stride, psd, rng)
    noff = i * stride
    ncpy = seglen if i == numseg - 1 else stride
    data.data.data[noff:noff+ncpy] = seg.data.data[0:ncpy]

# Injecting the generated noise on the waveform data
lalsim.SimAddInjectionREAL8TimeSeries(data, h, None)
FrWriteREAL8TimeSeries(data, 0)
frame_file_name = "X-{}-{}-{}-{}.gwf".format(channel, det, int(t0), int(seglen))
#fname = os.path.join(outdir, "{}-{}-{}-{}.gwf".format(detector, channel,int(t0), int(seglen)))
