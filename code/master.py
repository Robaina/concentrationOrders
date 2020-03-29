if __name__ == '__main__':

    from gurobipy import GRB
    import numpy as np
    import parameters as par
    import data
    import model

    GEM = data.prepareGEM(path_to_GEM=par.work_directory + par.model,
                          carbon_source=par.carbon_source,
                          uptake_rate=par.uptake_rate, biomass_threshold=par.beta,
                          loopless_fva=par.loopless_fva,
                          biomass_reaction_id=par.biomass_rxn_id)

    dG0_data = data.getFreeEnergyData(GEM,
                                      work_directory=par.work_directory,
                                      pH_i=par.pH_i,
                                      Ionic_strength=par.Ionic_strength,
                                      dG0_uncertainty_threshold=par.alpha)

    m, variables = model.buildMILPModel(GEM, dG0_data)
    # Fix internal pH (change to mM)
    logx_h_c = variables['logx']['logx_h_c']
    logx_h_c.lb = np.log(10**(-(par.pH_i + par.delta_pH) + 3))
    logx_h_c.ub = np.log(10**(-(par.pH_i - par.delta_pH) + 3))

    # Upper bound to internal o2 (mM)
    logx_o2_c = variables['logx']['logx_o2_c']
    logx_o2_c.ub = np.log(par.maxo2) # mM

    GEM_mets = [met.id for met in GEM.metabolites]
    met_p, met_q = GEM_mets[0], GEM_mets[1]

    obj_str = f'variables["logx"]["logx_{met_p}"] - variables["logx"]["logx_{met_q}"]'
    m.setObjective(eval(obj_str), GRB.MINIMIZE)
    m.setParam('OutputFlag', True)
    m.update()

    #m.optimize()
    # if m.status == GRB.Status.OPTIMAL:
    #     z = m.objval
    #     min_ratio = np.e**z
    #     print(min_ratio, m.SolCount)
    # else:
    #     print("fail")

    candidatePairs = model.findCandidatePairs(m, n_samples=100)
    print(len(candidatePairs))
    if len(candidatePairs) < 100:
        print(candidatePairs)
