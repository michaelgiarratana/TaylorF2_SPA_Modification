import numpy as np
import itertools
import bilby
import lalsimulation as lalsim
import lal
import argparse
import os
import glob
import sys
from astropy.cosmology import z_at_value, Planck18
from astropy import units
from gwpy.timeseries import TimeSeries
from lalframe import FrWriteREAL8TimeSeries

# Argument parsing setup
parser = argparse.ArgumentParser(description="Process parameter configuration file.")
parser.add_argument("--mass1", type=float, required=True, help="Mass of the primary source (in solar masses).")
parser.add_argument("--mass2", type=float, required=True, help="Mass of the secondary source (in solar masses).")
parser.add_argument("--s1", type=float, required=True, help="primary spin magnitude.")
parser.add_argument("--s2", type=float, required=True, help="secondary spin magnitude.")
parser.add_argument("--distance", type=float, required=True, help="Luminosity distance (in Mpc).")
parser.add_argument("--inclination", type=float, required=True, help="Injected inclination")
parser.add_argument("--eventT", type=float, required=True, help="Event time.")
parser.add_argument("--H_path", type=str, required=True, help="Path to H1 data.")
parser.add_argument("--L_path", type=str, required=True, help="Path to L1 data.")
parser.add_argument("--H_channel", type=str, required=True, help="Channel name for H1 data.")
parser.add_argument("--L_channel", type=str, required=True, help="Channel name for L1 data.")
parser.add_argument("--out_dir", type=str, required=True, help="Directory to hold output files.")
parser.add_argument("--run_label", type=str, required=True, help="Label for output files.")
parser.add_argument("--source", type=str, required=True, help="Type of binary (BBH, BNS, BHNS)")
args = parser.parse_args()  # Parse the arguments

# Derive other parameters from the input args
m1 = args.mass1
m2 = args.mass2
M = m1+m2
eta = (m1*m2)/(M**2)
a1 = args.s1
a2 = args.s2
dL = args.distance
chi_eff = (m1*a1 + m2*a2)/(M)


# Read the times series data from the gwf files
data_H = TimeSeries.read(args.H_path, channel = args.H_channel)
data_L = TimeSeries.read(args.L_path, channel = args.L_channel)

sampling_frequency = 4096
duration = len(data_H.data)/4096
outdir = args.out_dir
label = args.run_label

injection_parameters = dict(
    mass_1 = m1, # Msol
    mass_2 = m2,
    total_mass = M,
    symmetric_mass_ratio = eta,
    a_1 = a1,
    a_2 = a2,
    tilt_1=0.0, # Angles b/w the spin vectors of each black hole and the orbital angular momentum in radians
    tilt_2=0.0,
    phi_12=0.0, # The azimuthal angle b/w the two spin vectors in the plane perpendicular to the orbital angular momentum
    phi_jl=0, # Angle b/w the total angular momentum of the binary and the orbital angular momentum
    luminosity_distance = dL, # Mpc
    theta_jn = args.inclination,# Inclination angle b/w the line of sight and the orbital angular momentum of the binary, in radians
    phase=0, # The phase of the gravitational wave at the reference frequency
    ra=0, # The celestial coordinate defining the source's location on the sky, measured in radians
    dec=0, # The celestial coordinate defining the source's location on the sky, measured in radians
    geocent_time = args.eventT, #The GPS time at which the gravitational wave passes through the center of the Earth
    psi=0,# The polarization angle of the gravitational wave, describing the orientation of the wave’s polarization ellipse in the sky
)

# specify waveform arguments
waveform_arguments = dict(
    waveform_approximant="TaylorF2",  # waveform approximant name
    reference_frequency=0.0,  # Reference GW frequency (Hz). If 0Hz, reference point is coalescence
)

# Map string identifiers to the actual source models
source_model_map = {
    "BBH": bilby.gw.source.lal_binary_black_hole,
    "BNS": bilby.gw.source.lal_binary_neutron_star,
    "BHNS": bilby.gw.source.lal_binary_neutron_star 
}

# Set the source model based on gw source from args.source_type
frequency_domain_source_model = source_model_map.get(args.source_type.upper())

# Set up the waveform generator
waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
    sampling_frequency=sampling_frequency,
    duration=duration,
    frequency_domain_source_model=frequency_domain_source_model,
    parameters=injection_parameters,
    waveform_arguments=waveform_arguments,
)

# Combine the H1 and L1 data into a list to iterate over
all_data = [data_H, data_L]

# Set strain data for each detector
interferometers = bilby.gw.detector.InterferometerList(['H1', 'L1'])
for ifo, data in itertools.zip_longest(interferometers, all_data):
        ifo.strain_data.set_from_gwpy_timeseries(data)

# Create empty frequency series to accept psd data
psd = lal.CreateREAL8FrequencySeries(
    name = 'psd',
    epoch = 0,
    f0 = 15,
    deltaF = 1.0 / duration,
    sampleUnits = lal.SecondUnit,
    length = int(duration * sampling_frequency) // 2 + 1
)

# Create a frequency array from f_min (f0) to f_final - the nyquist frequency 
freqs = psd.f0 + np.arange(psd.data.length) * psd.deltaF
lalsim.SimNoisePSDaLIGOZeroDetHighPowerGWINC(psd, 15) # 15 is f_min in gwf construction

# Wrap into Bilby's PowerSpectralDensity format and assign to each interferometer
for ifo in interferometers:
    ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=freqs,
        psd_array=psd.data.data
    )

# Select the appropriate prior configuration based on source
source = args.source_type.upper()
priors = None  # Initialize priors

if source == "BBH":
    # Black Hole Binary prior setup
    # Set up all priors to be equal to a delta function at their values designated in the injection parameters
    priors = bilby.gw.prior.BBHPriorDict(injection_parameters.copy())

    # Replace fixed component masses with total_mass and eta
    priors.pop("mass_1")
    priors.pop("mass_2")

    priors['total_mass'] = bilby.core.prior.Uniform(name='total_mass', minimum=M * 0.5, maximum=M * 1.5)
    priors['symmetric_mass_ratio'] = bilby.core.prior.Uniform(name='symmetric_mass_ratio', minimum=0.1, maximum=0.25)
    priors['luminosity_distance'] = bilby.core.prior.PowerLaw(alpha=2, name='luminosity_distance',
                                                               minimum=dL * 0.5, maximum=dL * 1.5,
                                                               unit='Mpc', latex_label='$d_L$')
    priors['a_1'] = bilby.core.prior.Uniform(0, 0.99, 'a_1')
    priors['a_2'] = bilby.core.prior.Uniform(0, 0.99, 'a_2')

elif source == "BNS":
    # Neutron Star Binary prior setup
    # Set up all priors to be equal to a delta function at their values designated in the injection parameters
    priors = bilby.gw.prior.BNSPriorDict(injection_parameters.copy())

    # Replace fixed component masses with chirp mass and eta
    priors.pop("mass_1")
    priors.pop("mass_2")

    priors["chirp_mass"] = bilby.core.prior.Uniform(mc * 0.95, mc * 1.05, "chirp_mass")
    priors["mass_ratio"] = bilby.core.prior.Uniform(0.125, 1, "mass_ratio")
    priors['chi_1'] = bilby.core.prior.Uniform(-0.05, 0.05, 'chi_1')
    priors['chi_2'] = bilby.core.prior.Uniform(-0.05, 0.05, 'chi_2')

elif source == "BHNS":
    # Black Hole/Neutron Star binary prior setup (currently same as BNS)
    # Set up all priors to be equal to a delta function at their values designated in the injection parameters
    priors = bilby.gw.prior.BNSPriorDict(injection_parameters.copy())

    # Replace fixed component masses with chirp mass and eta
    priors.pop("mass_1")
    priors.pop("mass_2")

    priors["chirp_mass"] = bilby.core.prior.Uniform(mc * 0.95, mc * 1.05, "chirp_mass")
    priors["mass_ratio"] = bilby.core.prior.Uniform(0.125, 1, "mass_ratio")
    priors['chi_1'] = bilby.core.prior.Uniform(-0.05, 0.05, 'chi_1')
    priors['chi_2'] = bilby.core.prior.Uniform(-0.05, 0.05, 'chi_2')

else:
    raise ValueError(f"Unknown source type '{args.source_type}'. Please choose from 'BBH', 'BNS', or 'BHNS'.")

# Compute the likelihoods
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=interferometers, waveform_generator=waveform_generator
)

# Map source types to their corresponding conversion functions
conversion_function_map = {
    "BBH": bilby.gw.conversion.generate_all_bbh_parameters,
    "BNS": bilby.gw.conversion.generate_all_bns_parameters,
    "BHNS": bilby.gw.conversion.generate_all_bns_parameters,
}

# Select the appropriate conversion function based on the source
conversion_function = conversion_function_map.get(args.source_type.upper())


result = bilby.core.sampler.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler="dynesty",
    npoints=2000,
    walks = 100,
    maxmcmc = 5000,
    npool = 50,
    nact=10,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
    conversion_function=conversion_function,
    sample="unif",
    bootstrap=0,
    resume=True,
)

result.plot_corner()
