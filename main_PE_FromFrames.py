### Bilby script shared to me by Anarya Ray ###

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


def parse_cmd():
    parser=argparse.ArgumentParser()
    parser.add_argument("--luminosityDistance", help = "luminosity distance of injected source")
    parser.add_argument("--mass1", help = "mass of the primary component injected source")
    parser.add_argument("--mass2", help = "mass of the secondary component injected source")
    parser.add_argument("--waveform-approximant", help = "waveform approximant")
    parser.add_argument("--search-waveform-approximant", help = "search waveform approximant")
    parser.add_argument("--eos", help = "eos")
    parser.add_argument("--psd-dir", help = "directory containing power spectral densities of detectors")
    parser.add_argument("--ra",help= "injected right ascension")
    parser.add_argument("--dec",help= "injected declination")
    parser.add_argument("--inc",help= "injected inclination")
    parser.add_argument("--psi",help= "injected polarization")
    parser.add_argument("--observing-run", help='observing run')
    parser.add_argument("--chi1",help= "injected spin of primary componenet")
    parser.add_argument("--chi2",help= "injected spin of secondary componenet")
    parser.add_argument("--envName",help= "name of environment to change outdir path appropriately")
    args=parser.parse_args()
    return args
args=parse_cmd()

D=args.luminosityDistance
m1=float(args.mass1)
m2=float(args.mass2)
ra=float(args.ra)
dec=float(args.dec)
inc=float(args.inc)
psi=float(args.psi)
chi1 = float(args.chi1)
chi2 = float(args.chi2)
print(chi1,chi2)

envName = str(args.envName)

z=z_at_value(Planck18.luminosity_distance,float(D)*units.Mpc).value
obr=args.observing_run
approximant=args.waveform_approximant
search_approximant=args.search_waveform_approximant
eos=args.eos
psd_files=glob.glob(args.psd_dir+'*.txt')
mc=(m1*m2)**(3./5.)/(m1+m2)**(1./5.)#>1.9 or (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
mc*=(1.+z)
print(mc)
#sys.exit()

if mc<0.6 or mc>4:
    
    sys.exit()

#outdir = '/home/giarratana/outdir2'+approximant+'/'+eos+'/'+str(D)+'_'+str(int(m1*100)/100.)+'_'+str(int(m2*100)/100.)+'/'
outdir = '/home/giarratana/results/'
if(not os.path.exists(outdir)):
    os.makedirs(outdir)

file_to_det={'H1':"aligo",'L1':'aligo','V1':'avirgo','K1':'kagra'}
duration = 0

duration = 512

def lambda_from_piecewise(p,m):
    EOS=lalsim.SimNeutronStarEOS4ParameterPiecewisePolytrope(p[0],p[1],p[2],p[3])
    FAM=lalsim.CreateSimNeutronStarFamily(EOS)                                                                                                                                                                                         
    m_max=lalsim.SimNeutronStarMaximumMass(FAM)/lal.MSUN_SI
    m_max = int(m_max*1000)/1000
    rr = lalsim.SimNeutronStarRadius(m*lal.MSUN_SI, FAM)
    kk = lalsim.SimNeutronStarLoveNumberK2(m*lal.MSUN_SI, FAM)
    cc = m*lal.MRSUN_SI/rr
    return (2/3)*kk/(cc**5)

def lambda_from_spectral(p,m):
    EOS=lalsim.SimNeutronStarEOS4ParameterSpectralDecomposition(p[0],p[1],p[2],p[3])
    FAM=lalsim.CreateSimNeutronStarFamily(EOS)
    m_max=lalsim.SimNeutronStarMaximumMass(FAM)/lal.MSUN_SI
    m_max = int(m_max*1000)/1000
    rr = lalsim.SimNeutronStarRadius(m*lal.MSUN_SI, FAM)
    kk = lalsim.SimNeutronStarLoveNumberK2(m*lal.MSUN_SI, FAM)
    cc = m*lal.MRSUN_SI/rr
    return (2/3)*kk/(cc**5)


def lambda_from_name(m,name):
    EOS=lalsim.SimNeutronStarEOSByName(name)
    FAM=lalsim.CreateSimNeutronStarFamily(EOS)
    m_max=lalsim.SimNeutronStarMaximumMass(FAM)/lal.MSUN_SI
    m_max = int(m_max*1000)/1000
    rr = lalsim.SimNeutronStarRadius(m*lal.MSUN_SI, FAM)
    kk = lalsim.SimNeutronStarLoveNumberK2(m*lal.MSUN_SI, FAM)
    cc = m*lal.MRSUN_SI/rr
    return (2/3)*kk/(cc**5)

spectral_inj=np.array([+0.8651, +0.1548, -0.0151, -0.0002])
picewise_inj=np.array([])

if eos=='spectral':
    l1=lambda_from_spectral(spectral_inj,m1)
    l2=lambda_from_spectral(spectral_inj,m2)
elif eos=='piecewise':
    l1=lambda_from_piecewise(piecewise_inj,m1)
    l2=lambda_from_spectral(piecewise_inj,m2)
else:
    l1=lambda_from_name(m1,eos)
    l2=lambda_from_name(m2,eos)
print(bilby.gw.conversion.lambda_1_lambda_2_to_lambda_tilde(l1,l2,m1,m2))



#outdir = 'pe_dir'
label = 'bns_example'
bilby.core.utils.setup_logger(outdir=outdir, label=label)

minimum_frequency = 20
reference_frequency = 100
sampling_frequency = 4096.

data_H = TimeSeries.read("/home/giarratana/TaylorF2_SPA_Modification/X-Test400_H1-1197008864-16.gwf", channel="Test400-H1")
data_L = TimeSeries.read("/home/giarratana/TaylorF2_SPA_Modification/X-Test400_L1-1197008864-16.gwf", channel="Test400-L1")

all_data = [data_H, data_L]

# set up simulated data
np.random.seed(88170235)
injection_parameters = dict(
    chirp_mass=mc, mass_ratio=m2/m1, chi_1=chi1, chi_2=chi2, lambda_1=l1, lambda_2=l2,ra=ra, dec=dec, luminosity_distance=float(D), theta_jn=inc, psi=psi, phase=1.3, geocent_time=1187008882)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=dict(
        waveform_approximant=approximant,
        reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency
    )
)

if obr is None:
    obr='O4'

if obr=='O4':
    interferometers =bilby.gw.detector.InterferometerList(['H1', 'L1'])

    for ifo in interferometers:
        for fn in psd_files:
            if file_to_det[ifo.name] in fn:
                print(ifo.name,fn)
                ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(asd_file=fn)
elif obr=='O3':
    interferometers = bilby.gw.detector.InterferometerList(['H1', 'L1'])
    for ifo, data in itertools.zip_longest(interferometers, all_data):
        ifo.strain_data.set_from_gwpy_timeseries(data)



elif obr=='O2':
    interferometers =bilby.gw.detector.InterferometerList(['H1', 'L1'])
    for ifo, data in itertools.zip_longest(interferometers, all_data):
        ifo.strain_data.set_from_gwpy_timeseries(data)

else:
    print(obr)
    raise
#for interferometer in interferometers:
#    interferometer.minimum_frequency = minimum_frequency
#interferometers.set_strain_data_from_zero_noise(sampling_frequency, duration, start_time=ifo.start_time - duration + 2.)

#interferometers.inject_signal(
#    parameters=injection_parameters,
#    waveform_generator=waveform_generator
#)


# set up priors
priors = bilby.gw.prior.BNSPriorDict()
priors.pop("mass_1")
priors.pop("mass_2")
priors.pop("lambda_1")
priors.pop("lambda_2")

mc=injection_parameters['chirp_mass']
priors['chirp_mass'].minimum = mc * 0.95
priors['chirp_mass'].maximum = mc * 1.05
priors["mass_ratio"].minimum = 0.125
priors["mass_ratio"].maximum = 1
priors['chi_1'] = bilby.core.prior.Uniform(minimum=-0.05, maximum=0.05, name="chi_1")
priors['chi_2'] = bilby.core.prior.Uniform(minimum=-0.05, maximum=0.05, name="chi_2")
priors['lambda_tilde'] = bilby.core.prior.Uniform(0, 5000, name='lambda_tilde')
priors['delta_lambda'] = bilby.core.prior.Uniform(-5000, 5000, name='delta_lambda')
priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=injection_parameters['geocent_time'] - 0.1,
    maximum=injection_parameters['geocent_time'] + 0.1,
    name='geocent_time', latex_label='$t_c$', unit='$s$'
)
priors["luminosity_distance"] = bilby.core.prior.PowerLaw(alpha=2, name='luminosity_distance', minimum=min(10,float(D)-10), maximum=max(100,float(D)+100), unit='Mpc', latex_label='$d_L$')

# set up ROQ likelihood
search_waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby.gw.source.binary_neutron_star_roq,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
    waveform_arguments=dict(
        waveform_approximant=search_approximant,
        reference_frequency=reference_frequency
    )
)
roq_params = np.array(
    [(minimum_frequency, sampling_frequency / 2, duration, 0.8, 1.9, 0)],
    dtype=[("flow", float), ("fhigh", float), ("seglen", float), ("chirpmassmin", float), ("chirpmassmax", float), ("compmin", float)]
)
likelihood = bilby.gw.likelihood.ROQGravitationalWaveTransient(
    interferometers, search_waveform_generator, priors,
    linear_matrix="/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/simulations/roq/usedfor_uniformP_LTs/basis.hdf5", quadratic_matrix="/home/michael/projects/eos/GWXtreme_Tasks/year2/bilby_runs/simulations/roq/usedfor_uniformP_LTs/basis.hdf5", roq_params=roq_params,
    distance_marginalization=True, phase_marginalization=True
)
main_PE_FromFrames
# sampling
npool = 50
nact = 10
nlive = 2000
result = bilby.run_sampler(
    likelihood=likelihood, priors=priors, sampler='dynesty', use_ratio=True,
    nlive=nlive, walks=100, maxmcmc=5000, nact=nact, npool=npool,
    injection_parameters=injection_parameters, outdir=outdir, label=label, 
    conversion_function=bilby.gw.conversion.generate_all_bns_parameters,
    result_class=bilby.gw.result.CBCResult,
)    #


result.plot_corner()
