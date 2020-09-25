# Parameters
model = 'iML1515.json'
biomass_rxn_id = 'BIOMASS_Ec_iML1515_core_75p37M' # 'BIOMASS_Ec_iJO1366_core_53p95M'
carbon_source = 'EX_glc__D_e'
uptake_rate = 10  # mmol.min^-1.gDW^-1
biomass_threshold = beta = 0
x_min, x_max = 1e-8, 1e9  # mM
K_min, K_max = 1e-6, 100
gurobiFlag = True


# Output parameters
# fileName = 'graphData.json'
# dirName = 'cond6'
# directory = 'data/conditions/' + dirName
# notes = 'Growth on Glucose, no maximum metabolite sum constraint, including all pathways, reduced uncertainty in dG0 values, 0% of original'
work_directory = 'C:/Users/tinta/OneDrive/Documents/Projects/Ordering_of_Concentrations'
