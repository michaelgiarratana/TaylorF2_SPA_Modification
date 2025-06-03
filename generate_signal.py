import argparse
import sys
import math
import numpy as np
import lalsimulation as lalsim
import lal
from lalframe import FrWriteREAL8TimeSeries


# Argument parsing setup
parser = argparse.ArgumentParser(description="Process parameter configuration file.")
#parser.add_argument("--config", required=True, help="Path to parameter configuration file.")
parser.add_argument("--mass1", type=float, required=True, help="Mass of the primary source (in solar masses).")
parser.add_argument("--mass2", type=float, required=True, help="Mass of the secondary source (in solar masses).")
parser.add_argument("--S1x", type=float, required=True, help="x-component of the primary spin vector.")
parser.add_argument("--S1y", type=float, required=True, help="y-component of the primary spin vector.")
parser.add_argument("--S1z", type=float, required=True, help="z-component of the primary spin vector.")
parser.add_argument("--S2x", type=float, required=True, help="x-component of the secondary spin vector.")
parser.add_argument("--S2y", type=float, required=True, help="y-component of the secondary spin vector.")
parser.add_argument("--S2z", type=float, required=True, help="z-component of the secondary spin vector.")
parser.add_argument("--distance", type=float, required=True, help="Luminosity distance (in Mpc).")
parser.add_argument("--inclination", type=float, required=True, help="Injected inclination")
parser.add_argument("--phi_ref", type=float, required=True, help="Reference angle (in rad).")
parser.add_argument("--longAscNodes", type=float, required=True, help="Longitude of the ascending node of the binary's orbit.")
parser.add_argument("--meanPerAno", type=float, required=True, help="Mean anomaly at reference time.")
parser.add_argument("--eccentricity", type=float, required=True, help="Eccentricity of the binary's orbit.")
parser.add_argument("--deltaT", type=float, required=True, help="Time step (inverse of sampling frequency).")
parser.add_argument("--f_min", type=float, required=True, help="Starting frequency for waveform generation (Hz).")
parser.add_argument("--f_ref", type=float, required=True, help="Reference GW frequency (Hz).")
parser.add_argument("--LALparams", type=str, required=True, help="LAL dictionary used to pass additional parameters to the waveform generator.")
parser.add_argument("--approximant", required=True, help=" The waveform model.")
#
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

# Initialize parameter dictionary
params = {
    "m1": args.mass1 * lal.MSUN_SI,
    "m2": args.mass2 * lal.MSUN_SI,
    "S1x": args.S1x,
    "S1y": args.S1y,
    "S1z": args.S1z,
    "S2x": args.S2x,
    "S2y": args.S2y,
    "S2z": args.S2z,
    "distance": args.distance * 1e6 * lal.PC_SI,  # convert Mpc to meters
    "inclination": args.inclination,
    "phiRef": args.phi_ref,
    "longAscNodes": args.longAscNodes,
    "eccentricity": args.eccentricity,
    "meanPerAno": args.meanPerAno,
    "deltaT": float(args.deltaT),
    "f_min": args.f_min,
    "f_ref": args.f_ref,
    "LALparams": eval(args.LALparams),
    "approximant": lalsim.GetApproximantFromString(args.approximant),
}


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
#seglen = math.ceil(h.deltaT * h.data.length + 2*pad) # Computes the event duration in seconds (w/ a buffer)

# Below is the procedure for injecting noise from a psd onto the hp, hc data
segdur = np.ceil(len(h.data.data)/srate)*2 # segment duration
seglen = int(segdur * srate)
stride = seglen // 2  # stride between segments

# Below places the signal at the end of the noise
#start = segdur * np.floor(float(hp.epoch) / segdur)
#end = segdur * np.ceil((float(hp.epoch) + hp.data.length * hp.deltaT) / segdur)

# Below places the signal at the middle of the noise
#time_buffer = hp.epoch - (segdur * np.floor(float(hp.epoch) / segdur))
#start = hp.epoch - time_buffer/2
#end = segdur * np.ceil((float(hp.epoch) + hp.data.length * hp.deltaT) / segdur) + time_buffer/2

# Places the signal in the center of the noise NEED TO FIX FOR LONGER WAVEFORMS
start = eventT - segdur / 2
end = eventT + segdur / 2

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

#Adding the low pass filter and resampling the data
lal.LowPassREAL8TimeSeries(data,4096,0.9,8)

# Resample the data: downsample by a factor of 4 (original srate = 16384 -> 4096 Hz)
resampled_data = data.data.data[::4]

# Create a new REAL8TimeSeries for the resampled data
resampled_ts = lal.CreateREAL8TimeSeries(
    channel, data.epoch, data.f0, data.deltaT * 4, lal.StrainUnit, len(resampled_data)
    )
resampled_ts.data.data = resampled_data

FrWriteREAL8TimeSeries(resampled_ts, 0)
