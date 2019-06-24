# Constants (dG0 in kJ/mol) original dG0 data units
R = 8.3144598 * 1e-6 # kJ. K^-1. mmol^-1
T = 310.15 # 298.15 # K

# Parameters
number_preprocessing_samples = 10
loopless_fva = False
numerator_metabolites = p_mets = None
denominator_metabolites = q_mets = None
dG0_uncertainty_threshold = alpha = 1 # filter noisy dG0 data, percentage of mean
carbon_source = 'EX_glc__D_e'
uptake_rate = 10 # mmol.min^-1.gDW^-1
biomass_threshold = beta = 0.95
dG0_error_fraction = gamma = 1 # allowed dG0 deviation from measured values: +- gamma * measured dG0
Gibbs_eps = 1e-8 # kJ/mmol
M = 1e8 # big number to unconstrain reactions with no dG0 data
x_min, x_max = 1e-4, 2e2 # mM
Intracellular_pH = pH_i = 7.4
pH_deviation = delta_pH = 0.3
Ionic_strength = 0.25
maximum_internal_O2 = maxo2 = 0.25 # mM  {uM = 250 (Potzkei et al, 2012)}
# Fluxes in mmol.gDW^-1.min^-1
# (in Teppe et al 2013* they use 1e-5, 1e-1!
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)
fileName = 'graphData.json'
dirName = 'Cond11'
directory ='data/conditions/' + dirName
notes = 'Growth on Glucose, constrained maximum [o2]'
work_directory = 'C:/Users/tinta/OneDrive/Documents/Projects/Ordering_of_Concentrations/'
