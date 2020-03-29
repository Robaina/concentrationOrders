# Constants (dG0 in kJ/mol) original dG0 data units
R = 8.3144598  # kJ. K^-1. mmol^-1
T = 310.15  # 298.15 # K

# Parameters
model = 'iML1515.json'
# model = 'iJO1366.json'
biomass_rxn_id = 'BIOMASS_Ec_iML1515_core_75p37M'
# biomass_rxn_id = 'BIOMASS_Ec_iJO1366_core_53p95M'
number_preprocessing_samples = 1
numerator_metabolites = p_mets = None
denominator_metabolites = q_mets = None
dG0_uncertainty_threshold = alpha = 1  # filter noisy dG0 data, percentage of mean
carbon_source = 'EX_glc__D_e'
uptake_rate = 10  # mmol.min^-1.gDW^-1
biomass_threshold = beta = 1
dG0_error_fraction = gamma = 0.1  # allowed dG0 deviation from measured values: +- gamma * measured dG0
dG0_eps = 10  # kJ/mmol
# v_eps = 1e-3 # mmol.min^-1.gDW^-1
M = 1e8  # big number to unconstrain reactions with no dG0 data
x_min, x_max = 1e-3, 1e1  # mM
max_sum_x = 350 # mM
min_sum_x = 300 # mM
# min_sum_logx = -7550 # this is minimum of sum_i log(x_i)
# max_sum_logx = -7550
Intracellular_pH = pH_i = 7.4
pH_deviation = delta_pH = 0.4
Ionic_strength = '200 mM'
maximum_internal_O2 = maxo2 = 0.25  # mM  {uM = 250 (Potzkei et al, 2012)}
# Real-time determination of intracellular oxygen in bacteria using a genetically encoded FRET-based biosensor
# Fluxes in mmol.gDW^-1.min^-1
# (in Teppe et al 2013* they use 1e-5, 1e-1
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)


fileName = 'graphData.json'
dirName = 'Cond11'
directory = 'data/conditions/' + dirName
notes = 'Growth on Glucose, constrained maximum [o2]'
work_directory = 'C:/Users/tinta/OneDrive/Documents/Projects/Ordering_of_Concentrations'
