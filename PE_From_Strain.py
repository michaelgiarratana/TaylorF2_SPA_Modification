# PE from standard method of injecting the parameters
import numpy as np
import bilby
import lalsimulation as lalsim
import lal
import os
import glob
import sys
from astropy.cosmology import z_at_value, Planck18
from astropy import units
from gwpy.timeseries import TimeSeries
from lalframe import FrWriteREAL8TimeSeries
import pickle
import argparse

# New PE code where bilby injects the modified hp, hc data from a pickle file.

# Argument parsing setup
parser = argparse.ArgumentParser(description="Process parameter configuration file.")
parser.add_argument("--mass1", type=float, required=True, help="Mass of the primary source (in solar masses).")
parser.add_argument("--mass2", type=float, required=True, help="Mass of the secondary source (in solar masses).")
parser.add_argument("--s1", type=float, required=True, help="primary spin magnitude.")
parser.add_argument("--s2", type=float, required=True, help="secondary spin magnitude.")
parser.add_argument("--distance", type=float, required=True, help="Luminosity distance (in Mpc).")
parser.add_argument("--inclination", type=float, required=True, help="Injected inclination")
parser.add_argument("--eventT", type=float, required=True, help="Event time.")
parser.add_argument("--phi_ref", type=float, required=True, help="Reference GW phase (radians)")
parser.add_argument("--f_min", type=float, required=True, help="Starting frequency for waveform generation (Hz).")
parser.add_argument("--f_max", type=float, required=True, help="Ending frequency for waveform generation (Hz).")
parser.add_argument("--f_ref", type=float, required=True, help="Reference GW frequency (Hz).")
parser.add_argument("--out_dir", type=str, required=True, help="Directory to hold output files.")
parser.add_argument("--run_label", type=str, required=True, help="Label for output files.")
parser.add_argument("--source_type", type=str, required=True, help="Type of binary (BBH, BNS, BHNS)")
parser.add_argument("--injection_type", type=str, required=True, help="Method of signal injection (bilby or pickle)")
parser.add_argument("--strain_path", type=str, required=True, help="Path to strain pickle file.")

args = parser.parse_args()  # Parse the arguments


# Set the duration and sampling frequency of the data segment that we're
# going to inject the signal into

if args.source_type == 'BBH':
    duration = 32
elif args.source_type == "BNS":
    duration = 128*2

#duration = 32.0 # for BBH
#duration = 128.0 # for BNS
sampling_frequency = 4096.0

# Specify the output directory and the name of the simulation.
outdir = "/home/giarratana/TaylorF2_SPA_Modification/results/"
label = args.run_label
#bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Set up a random seed for result reproducibility
bilby.core.utils.random.seed(88170235)

# Pulling injection parameters from parser
m1 = args.mass1
m2 = args.mass2
q = m2/m1
M = m1+m2
eta = (m1*m2)/(M**2)
mc = ((m1*m2)**(3/5))/((m1+m2)**(1/5))
a1 = args.s1
a2 = args.s2
dL = args.distance
chi_eff = (m1*a1 + m2*a2)/(M)

# Defining injection params dict
injection_parameters = dict(
    mass_1 = m1,
    mass_2 = m2,
    chirp_mass = mc,
    mass_ratio = q,
    total_mass = M,
    symmetric_mass_ratio = eta,
    a_1 = a1,
    a_2 = a2,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    luminosity_distance = dL,
    theta_jn=0.0,
    psi=0,
    phase = args.phi_ref,
    geocent_time = args.eventT,
    ra=0,
    dec=0,
)

# Slight modification to the injections params for BNS
if args.source_type == 'BNS':
    # Replace a_1, a_2 with chi_1, chi_2
    injection_parameters['chi_1'] = injection_parameters.pop('a_1')
    injection_parameters['chi_2'] = injection_parameters.pop('a_2')

    # Add tidal deformabilities
    injection_parameters['lambda_1'] = 0
    injection_parameters['lambda_2'] = 0


# Arguments passed into the source model
waveform_arguments = dict(
    waveform_approximant = "TaylorF2",
    reference_frequency = args.f_ref,
    minimum_frequency = args.f_min,
    maximum_frequency = args.f_max,
)

# Map string identifiers to the actual source models
source_model_map = {
        "BBH": bilby.gw.source.lal_binary_black_hole,
        "BNS": bilby.gw.source.lal_binary_neutron_star,
        #"BNS": bilby.gw.source.binary_neutron_star_roq,
        "BHNS": bilby.gw.source.lal_binary_neutron_star
}

# Set the source model based on gw source from args.source_type
frequency_domain_source_model = source_model_map.get(args.source_type.upper())

# Create the waveform_generator using the specified  source function
waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=frequency_domain_source_model,
    #parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)

if args.source_type == "BBH":
    begin = 4, # for BBH
elif args.source_type == "BNS":
    begin = 100, # for BNS


# Set up interferometers and the psd/asd from a corresponding txt file
ifos = bilby.gw.detector.InterferometerList(["H1", "L1"])
ifos[0].power_spectral_density = bilby.gw.detector.PowerSpectralDensity(asd_file = "/home/giarratana/TaylorF2_SPA_Modification/aLIGO_ZERO_DET_high_P_asd.txt")
ifos[1].power_spectral_density = bilby.gw.detector.PowerSpectralDensity(asd_file = "/home/giarratana/TaylorF2_SPA_Modification/aLIGO_ZERO_DET_high_P_asd.txt")
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 4, 
)

# Define the signal injection method
injection = args.injection_type


# Set up the injection call based on the injection_type
if injection == "pickle":
    
    # Load the waveform_dict pickle file
    with open(args.strain_path, "rb") as f:
        waveform_dict=pickle.load(f)

    ifos.inject_signal(
            injection_polarizations=waveform_dict, 
            parameters=injection_parameters
            )
elif injection == "bilby":
    ifos.inject_signal(
            waveform_generator=waveform_generator, 
            parameters=injection_parameters
            )
else: 
    raise ValueError(f"Unknown injection type: {injection}") 




# Set up a PriorDict, which inherits from dict.
# By default we will sample all terms in the signal models.  However, this will
# take a long time for the calculation, so for this example we will set almost
# all of the priors to be equall to their injected values.  This implies the
# prior is a delta function at the true, injected value.  In reality, the
# sampler implementation is smart enough to not sample any parameter that has
# a delta-function prior.
# The above list does *not* include mass_1, mass_2, theta_jn and luminosity
# distance, which means those are the parameters that will be included in the
# sampler.  If we do nothing, then the default priors get used.

# Select the appropriate prior configuration based on source
source = args.source_type.upper()
priors = None  # Initialize priors

if source == 'BBH':
    priors = bilby.gw.prior.BBHPriorDict(injection_parameters.copy())
    priors['geocent_time'] = bilby.core.prior.Uniform(name='geocent_time', minimum=injection_parameters['geocent_time'] - 0.1, maximum=injection_parameters['geocent_time'] + 0.1)
    priors['phase'] = bilby.core.prior.Uniform(name='phase', minimum=0, maximum=np.pi)
    
    #priors['mass_1'] = bilby.core.prior.Uniform(name='mass_1', minimum=20, maximum=45)
    #priors['mass_2'] = bilby.core.prior.Uniform(name='mass_2', minimum=20, maximum=45)
    
    priors.pop("mass_1")
    priors.pop("mass_2")
    priors.pop("symmetric_mass_ratio")
    priors.pop("total_mass")
    priors['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=mc*0.95, maximum=mc*1.05)
    priors['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.1, maximum=1)
    
    priors.pop('luminosity_distance')
    priors['luminosity_distance'] = bilby.core.prior.PowerLaw(alpha=2, name='luminosity_distance',minimum=100, maximum=500, unit='Mpc', latex_label='$d_L$')
    #priors.pop("mass_1")
    #priors.pop("mass_2")
    #priors.pop("mass_ratio")
    #priors.pop("chirp_mass")
    #priors['total_mass'] = bilby.core.prior.Uniform(name='total_mass', minimum=8, maximum=20)
    #priors['symmetric_mass_ratio'] = bilby.core.prior.Uniform(name='symmetric_mass_ratio', minimum=0.1, maximum=0.7)
    #priors['a_1'] = bilby.core.prior.Uniform(0, 0.9, 'a_1')
    #priors['a_2'] = bilby.core.prior.Uniform(0, 0.9, 'a_2')

elif source == "BNS":
    # Neutron Star Binary prior setup
    # Set up all priors to be equal to a delta function at their values designated in the injection parameters
    priors = bilby.gw.prior.BNSPriorDict(injection_parameters.copy())
    priors['geocent_time'] = bilby.core.prior.Uniform(name='geocent_time', minimum=injection_parameters['geocent_time'] - 0.1, maximum=injection_parameters['geocent_time'] + 0.1)
    priors['phase'] = bilby.core.prior.Uniform(name='phase', minimum=0, maximum=np.pi)
    
    #priors['mass_1'] = bilby.core.prior.Uniform(name='mass_1', minimum=20, maximum=45)
    #priors['mass_2'] = bilby.core.prior.Uniform(name='mass_2', minimum=20, maximum=45)
    
    #priors.pop('luminosity_distance')
    #priors['luminosity_distance'] = bilby.core.prior.PowerLaw(alpha=2, name='luminosity_distance',minimum=20, maximum=100, unit='Mpc', latex_label='$d_L$')

    priors.pop("mass_1")
    priors.pop("mass_2")
    priors.pop("symmetric_mass_ratio")
    priors.pop("total_mass")
    priors['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=1.19, maximum=1.39)
    priors['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.3, maximum=1)
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


# Perform a check that the prior does not extend to a parameter space longer than the data
priors.validate_prior(duration, args.f_min)


# Compute the appropriate likelihoods based on the source type
if args.source_type.upper() in ["BNS", "BBH", "NSBH"]:
    # Initialise the likelihood by passing in the interferometer data (ifos) and the waveform generator
    # Standard likelihood for BBH
    likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=waveform_generator,
        priors=priors, #For the likelihood marginalization
        distance_marginalization=True,
       # time_marginalization=True
    )
elif args.source_type.upper() in []:#"BNS", "NSBH"]: # Currently not in use - returned poor recovery for unknown reason
    # ROQ params for BNS likelihood
    roq_params = np.array(
    [(args.f_min, sampling_frequency / 2, duration, 0.8, 1.9, 0)],
    dtype=[("flow", float), ("fhigh", float), ("seglen", float), ("chirpmassmin", float), ("chirpmassmax", float), ("compmin", float)]
    )
    # ROQ likelihood for BNS/NSBH
    likelihood = bilby.gw.likelihood.ROQGravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=waveform_generator,
        priors=priors,
        linear_matrix="/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/simulations/roq/usedfor_uniformP_LTs/basis.hdf5",
        quadratic_matrix="/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/simulations/roq/usedfor_uniformP_LTs/basis.hdf5",
        roq_params=roq_params,
        distance_marginalization=True,
        phase_marginalization=False
    )
else:
    raise ValueError(f"Unknown source type: {args.source_type}")

# Defining custom conversion functions with custom waveform defaults
from bilby.gw.conversion import _generate_all_cbc_parameters, convert_to_lal_binary_black_hole_parameters, convert_to_lal_binary_neutron_star_parameters, generate_tidal_parameters
def generate_all_bbh_parameters_custom(sample, likelihood = None, priors = None, npool = 1):
    waveform_defaults = {
            'reference_frequency': args.f_ref, 'waveform_approximant': 'TaylorF2',
            'minimum_frequency': args.f_min}
    output_sample = _generate_all_cbc_parameters(
            sample, defaults = waveform_defaults,
            base_conversion = convert_to_lal_binary_black_hole_parameters,
            likelihood = likelihood, priors = priors, npool = npool)
    return output_sample

def generate_all_bns_parameters_custom(sample, likelihood = None, priors = None, npool = 1):
    waveform_defaults ={
            'reference_frequency': args.f_ref, 'waveform_approximant': 'TaylorF2',
            'minimum_frequency': args.f_min}
    output_sample = _generate_all_cbc_parameters(
            sample, defaults = waveform_defaults,
            base_conversion = convert_to_lal_binary_neutron_star_parameters,
            likelihood = likelihood, priors = priors, npool = npool)
    try:
        output_sample = generate_tidal_parameters(output_sample)
    except KeyError as e:
        logger.debug(
                "Generation of tidal parameters failed with message {}".format(e))
        return output_sample

# Map source types to their corresponding conversion functions
#conversion_function_map = {
#    "BBH": bilby.gw.conversion.generate_all_bbh_parameters,
#    "BNS": bilby.gw.conversion.generate_all_bns_parameters,
#    "BHNS": bilby.gw.conversion.generate_all_bns_parameters,
#}

conversion_function_map = {
        "BBH": generate_all_bbh_parameters_custom,
        "BNS": generate_all_bns_parameters_custom,
        "BHNS": generate_all_bns_parameters_custom,
        }

# Select the appropriate conversion function based on the source
conversion_function = conversion_function_map.get(args.source_type.upper())

result = bilby.core.sampler.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler="dynesty",
    npoints=3000,
    walks = 100,
    maxmcmc = 5000,
    npool = 60,
    nact=10,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
    #conversion_function=conversion_function,
    sample="unif",
    bootstrap=0,
    resume=True,
)

result.plot_corner()
