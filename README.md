# Ordering of metabolite concentrations in _E. coli_
This is an ongoing project in which I use an extended version of thermodynamic flux balance analysis to find metabolite pairs, $i,j$, such that:
$$x_i \geq x_j$$
in any steady state of the system, where $x_i$ is the concentration of metabolite $i$.

I'm using data on $\Delta G$ from the equilibrator api.

---

  Several studies have shown that certain metabolites maintain an order relation in their concentrations under different steady states of \emph{Escherichia coli} \cite{Bennett2008,Bennett2009a}. This order relation may be the result of diverse constraints operating at steady state, such as stoichiometric, thermodynamic and growth constraints. In the following, we will derive a theoretical explanation for this ordering of concentrations.

  We describe the concentration dynamics of the biochemical network with the general system,

  \begin{equation}
    \label{eq:1}
    \frac{dx_i}{dt} = \sum_j n_{ij} v_j(x),
  \end{equation}

  \noindent where, $n_{ij}$ represents the stoichiometric coefficient of metabolite $i$ in reaction $j$, with $n_{ij} < 0$ if it is a substrate of the reaction, $n_{ij} > 0$ if it is a product, $x$ the metabolite concentrations and $v(x)$ the metabolic fluxes, which are a function of the metabolite concentration --- the exact form provided by the selected kinetic law. All reactions in \ref{eq:1} are reversible, with the exception of the biomass production (pseudo)reaction, and the forward and backward direction are represented as two different reactions. Thus, the net flux of a reaction $v_j = v_j^{for} - v_j^{back}$. We assume that cells are growing at steady state, thus the flux through the biomass reaction, $v_{bio} > \gamma v^{max}_{bio}$, with $\gamma \in [0, 1]$, a fraction of the theoretical maximum.

  Thus far, we have determined that certain reactions are irreversible in an escenario where cells grow at steady state. Additionally, the second law of thermodynamics imposes that flux of free energy $g_j = \Delta^{\circ} G_j v_j < 0$ for a reaction to have non-zero flux \cite{Kondepudi2014a}. In our case, we have already determined the direction of the reaction with the linear programs in \ref{eq:2}, hence we only need the reaction Gibbs free energy

  \begin{equation}
    \label{eq:2}
    \Delta^{\circ} G_j < 0,
  \end{equation}

  \noindent where,

  \begin{equation}
    \label{eq:3}
    \Delta^{\circ} G_j = \Delta^{\circ} G_{r(j)} + RT \sum_i n_{ij} \log{x_i}
  \end{equation}

  \noindent in which $\Delta^{\circ} G_{r(j)} = \sum_i n_{ij}\Delta^{\circ} G_{f(i)}$. Further, the energies of formation $\Delta^{\circ} G_{f(i)}$ of the metabolites participating in the reaction can be estimated with the component contribution method \cite{Noor2013}.


```python
from gurobipy import GRB
from importlib import reload
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import parameters as par
import data
import model
import output
reload(par)

import pickle

def saveToPickleFile(python_object, path_to_file='object.pkl'):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file,'wb')
    pickle.dump(python_object, out_file)
    out_file.close()
    
def readFromPickleFile(path_to_file='object.pkl'):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file,'rb')
    python_object = pickle.load(in_file)
    return python_object


GEM = data.prepareGEM(path_to_GEM=f'{par.work_directory}/{par.model}',
                      carbon_source=par.carbon_source,
                      uptake_rate=par.uptake_rate,
                      loopless=False,
                      biomass_reaction_id=par.biomass_rxn_id)

dG0_data = data.getFreeEnergyData(GEM,
                                  work_directory=par.work_directory,
                                  pH_i=par.pH_i,
                                  Ionic_strength=par.Ionic_strength,
                                  dG0_uncertainty_threshold=par.alpha)

GEM
```

    Downloading package metadata...
    Fragments already downloaded
    Downloading package metadata...
    Fragments already downloaded
    Using license file C:\Users\tinta\gurobi.lic
    Academic license - for non-commercial use only
    Read LP format model from file C:\Users\tinta\AppData\Local\Temp\tmpitt8bst0.lp
    Reading time = 0.05 seconds
    : 1147 rows, 4276 columns, 17168 nonzeros
    





<table>
    <tr>
        <td><strong>Name</strong></td>
        <td>iML1515</td>
    </tr><tr>
        <td><strong>Memory address</strong></td>
        <td>0x0274bd05c848</td>
    </tr><tr>
        <td><strong>Number of metabolites</strong></td>
        <td>1147</td>
    </tr><tr>
        <td><strong>Number of reactions</strong></td>
        <td>2138</td>
    </tr><tr>
        <td><strong>Number of groups</strong></td>
        <td>0</td>
    </tr><tr>
        <td><strong>Objective expression</strong></td>
        <td>1.0*BIOMASS_Ec_iML1515_core_75p37M - 1.0*BIOMASS_Ec_iML1515_core_75p37M_reverse_35685</td>
    </tr><tr>
        <td><strong>Compartments</strong></td>
        <td>cytosol, extracellular space, periplasm</td>
    </tr>
  </table>



## 2 iML1515 model and gibbs data with some examples


```python
# Let's look at the Cytidine kinase
rxn_id = 'CYTDK2'
rxn = GEM.reactions.get_by_id('CYTDK2')
rxn
```





<table>
    <tr>
        <td><strong>Reaction identifier</strong></td><td>CYTDK2</td>
    </tr><tr>
        <td><strong>Name</strong></td><td>Cytidine kinase (GTP)</td>
    </tr><tr>
        <td><strong>Memory address</strong></td>
        <td>0x022b10302f08</td>
    </tr><tr>
        <td><strong>Stoichiometry</strong></td>
        <td>
            <p style='text-align:right'>cytd_c + gtp_c --> cmp_c + gdp_c + h_c</p>
            <p style='text-align:right'>Cytidine + GTP C10H12N5O14P3 --> CMP C9H12N3O8P + GDP C10H12N5O11P2 + H+</p>
        </td>
    </tr><tr>
        <td><strong>GPR</strong></td><td>b2066</td>
    </tr><tr>
        <td><strong>Lower bound</strong></td><td>0.0</td>
    </tr><tr>
        <td><strong>Upper bound</strong></td><td>1000.0</td>
    </tr>
</table>





```python
# Its estimated Gibbs free energy of reaction is
dG0, error = dG0_data[rxn_id]["dG0"], dG0_data[rxn_id]["error"]
print(f'Reaction {rxn.name} with dG_0 = {dG0:.2f} +/- {error:.2f} kJ/mmol')
```

    Reaction Cytidine kinase (GTP) with dG_0 = -11.50 +/- 3.10 kJ/mmol
    

## 3 Constructing the MILP in Gurobi

  Now, the logarithm is a monotonically increasing function, hence $\log{x_p} > \log{x_q} \implies x_p > x_q$. Furthermore, we can establish if $\log{x_p} > \log{x_q}$ with the following convex optimization problem ($\mathrm{OP}_1$):
 
 
\begin{align}
    \begin{aligned}
      \label{eq:4}
      &v_{bio}^* = \max_{\substack{ \log{x} \in \rm I\!R^m, \\
                            v \in \rm I\!R^n_{\geq 0}, \\
                            \Delta^{\circ} G_{r} \in \rm I\!R^k, \\
                            y \in \{0, 1\}^k }} \; v_{bio}
      \\
      &\mathrm{s.t.}
      \\
      &1.\;Sv = 0
      \\
      &2.\;\Delta^{\circ} G_{r(j)} + RT \sum_i n_{ij} \log{x_i} - (1 - y_j^{(+, -)})M< 0 \; \forall j \in \mathrm{R_{\Delta^{\circ} G_{r}}}
      \\
      &3.\;\Delta^{\circ} G_{r}^{\dagger} - \epsilon \leq \Delta^{\circ} G_{r} \leq \Delta^{\circ} G_{r}^{\dagger} + \epsilon
      \\
      &4.\;\sum_i x_i \leq \; 300 \,\mathrm{mM}
      \\
      &5.\;\log{x}_{min} \leq \log{x} \leq \log{x}_{max}
      \\
      &6.\;v_{min} \leq v \leq v_{max}
      \\
      &7.\;v^{(+, -)}_j \leq y_j^{(+, -)} v_{max(j)} \; \forall j \in \mathrm{R_{\Delta^{\circ} G_{r}}}
      \\
      &8.\;y^+ + y^- \leq 1
    \end{aligned}
\end{align}

Here $\mathrm{OP}_2$


\begin{align}
\begin{aligned}
  \label{eq:5}
  &z = \min_{\substack{ \log{x} \in \rm I\!R^m, \\
                        v \in \rm I\!R^n_{\geq 0}, \\
                        \Delta^{\circ} G_{r} \in \rm I\!R^k, \\
                        y \in \{0, 1\}^k }} \; \log{x_p} - \log{x_q}
  \\
  &\mathrm{s.t.}
  \\
  &(1-8. \; \mathrm{OP}_1)
  \\
  &9.\;v_{bio} \geq \gamma v^*_{bio}
\end{aligned}
\end{align}

__NOTE__
What if Kms are already optimized to guarantee that metabolite concentrations that are thermodynamically feasible are already enzyme saturating? i.e., the system (or at least part of the system) has linear kinetics which are concentration independent. 


```python
# Build gurobi TFBA model
reload(par)
m, variables = model.buildMILPModel(GEM, dG0_data)

# Fix internal pH (change to mM)
logx_h_c = variables['logx']['logx_h_c']
logx_h_c.lb = np.log(10**(-(par.pH_i + par.delta_pH) + 3))
logx_h_c.ub = np.log(10**(-(par.pH_i - par.delta_pH) + 3))

# Upper bound to internal h2o (mM)
logx_h2o_c = variables['logx']['logx_h2o_e']
logx_h2o_c.ub = np.log(1e2) # mM

# Upper bound to internal o2 (mM)
logx_o2_c = variables['logx']['logx_o2_c']
logx_o2_c.ub = np.log(par.maxo2) # mM

# Fixed external glucose concentration
logx_glu_e = variables['logx']['logx_glu__L_e']
logx_glu_e.lb = np.log(20) # mM 22.2
logx_glu_e.ub = np.log(23) # mM

# Find biomass maximum under thermodynamic constraints and constraint biomass reaction
m.setObjective(variables['v']['v_' + par.biomass_rxn_id], GRB.MAXIMIZE)
m.update()
m.setParam('OutputFlag', False)
m.optimize()

try:
    print(f'Maximum growth rate: {m.objval:.2f} h^-1')
    max_biomass_flux = m.getVarByName('v_' + par.biomass_rxn_id).x
    m.getVarByName('v_' + par.biomass_rxn_id).lb = par.beta * max_biomass_flux
except:
    print('Model is infeasible!')
```

    Maximum growth rate: 0.88 h^-1
    


```python
# # Some sanity checks of the solution
# conc_sum = 0
# for var in m.getVars():
#     if 'x_' in var.varName[:2]:
# #         conc_sum += var.x
# total_conc_internal_mets = conc_sum - m.getVarByName('x_glu__L_e').x
# print(f'Total concentration of internal metabolites: {total_conc_internal_mets:.2f} mM')
# print(f'Total sum of metabolite concentrations: {conc_sum:.2f} mM')
```

## 4 Finding candidate ordered pairs

\begin{align}
\begin{aligned}
  \label{eq:10}
  &x^* = \mathrm{argmin}_{\substack{ \log{x} \in \rm I\!R^m, \\
                        v \in \rm I\!R^n_{\geq 0}, \\
                        \Delta^{\circ} G_{r} \in \rm I\!R^k, \\
                        y \in \{0, 1\}^k }} \; \sum_{i=1}^m {(\epsilon^+_i + \epsilon^-_i)}
  \\
  &\mathrm{s.t.}
  \\
  &(1-8. \; \mathrm{OP}_1)
  \\
  &(9. \; \mathrm{OP}_2)
  \\
  &10.\;\epsilon^+ - \epsilon^- = \log{x_{rand}} - \log{x}
  \\
  &11\;\epsilon^+, \epsilon^- \geq 0
\end{aligned}
\end{align}


```python
# Find candidate ordered metabolite pairs
import datetime
reload(par)
m.update()
start = datetime.datetime.now()
X_sample, candidatePairs, sampling_model = model.findCandidatePairs(m, n_samples=par.n_samples)
end = datetime.datetime.now()
print (f'Ellapsed time: {end - start}')
print(f'There are {len(candidatePairs)} candidate pairs')
```

    Ellapsed time: 0:03:30.377159
    There are 3647 candidate pairs
    

### 4.2 Restricting search to central metabolism

We can also look for pairs that do not involved constrained metabolites, like h_c, o_c and glu__L_c. I'm sure this will reduce the number of pairs to evaluate considerably.


```python
# Remove candidate pairs involving constrained metabolites
constrained_mets = ['h_c', 'o2_c', 'h2o_e', 'glu__L_e', 'h2o_c', 
                    'h_e', 'co2_c', 'co2_e', 'o2_e', 'h2o_p', 'o2_p', 
                    'co2_p', 'h_p', 'ppi_e', 'ppi_c', 'ppi_p', 
                    'pi_c', 'pi_e', 'pi_p']
filtered_pairs = [pair for pair in candidatePairs 
                  if (pair[0] not in constrained_mets 
                      and pair[1] not in constrained_mets)]

print(f'There are {len(filtered_pairs)} filtered candidate pairs')
```

    There are 287 filtered candidate pairs
    


```python
systems_df = pd.read_excel(f'{par.work_directory}/iML1515_subsystems.xlsx')
central_metabolism = ['Carbohydrate metabolism'
                      ,'Amino acid metabolism', 'Nucleotide metabolism',
                       'Lipid metabolism']
central_metabolites = data.findMetabolitesInPathways(GEM, systems_df,
                                                     central_metabolism)

centralCandidatePairs = []
for pair in filtered_pairs:
    met_i, met_j = pair
    if met_i in central_metabolites and met_j in central_metabolites:
        centralCandidatePairs.append(pair)

print(f'There are {len(centralCandidatePairs)} filtered and central candidate pairs')
```

    There are 163 filtered and central candidate pairs
    


```python
saveToPickleFile(X_sample, 'X_sample.pkl')
```

## 4.3 Evaluating candidate ordered pairs


```python
# Evaluate metabolite orders from candidate pairs
import datetime
start = datetime.datetime.now()
ordered_pairs = model.findConcentrationOrderedPairs(sampling_model, filtered_pairs)
end = datetime.datetime.now()
print (f'Ellapsed time: {end - start}')
```

    Evaluating concentration-ordered pairs...
    Ellapsed time: 0:06:29.622865
    


```python
print(f'There are {len(ordered_pairs)} true ordered pairs')
ordered_pairs
```

    There are 9 true ordered pairs
    




    [['dhap_c', 'g3p_c', 1.0021431070183684],
     ['3pg_c', '2pg_c', 1.001611461402214],
     ['5caiz_c', '5aizc_c', 1.012661974412411],
     ['mlthf_c', '5mthf_c', 98.64911545693387],
     ['mlthf_c', 'nad_c', 98.64911545693307],
     ['g6p_c', 'f6p_c', 1.0009747501109063],
     ['gam6p_c', 'gam1p_c', 1.0035629977225256],
     ['nadh_c', '5mthf_c', 98.64911545693369],
     ['nadh_c', 'nad_c', 98.64911545693307]]




```python
# Retrieve ordered pairs from pickle files
# import os
# import parameters as par
# directory=par.work_directory+'/'+par.directory+'/ordered_pairs'
# files = os.listdir(directory)
# ordered_pairs = []
# for file in files:
#     pairs = readFromPickleFile(directory+'/'+file)
#     ordered_pairs += pairs
```


```python
# ordered_pairs
```

## 4.4 Evaluating concentration order DAG


```python
reload(par)
reload(output)
A_plus, total_mets = output.buildDAGadjacencyMatrix(GEM, ordered_pairs)
reduced_met_orders = output.removeEdgesToProton(A_plus, total_mets)
graphData = output.buildGraphJSONData(GEM, reduced_met_orders)
output.writeOutputFiles(graphData)
```

     
    Writing files...
    


```python
from gurobipy import GRB
from importlib import reload
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import parameters as par
import data
import model
import output
reload(par)

GEM = data.prepareGEM(path_to_GEM=f'{par.work_directory}/{par.model}',
                      carbon_source=par.carbon_source,
                      uptake_rate=par.uptake_rate,
                      loopless=False,
                      biomass_reaction_id=par.biomass_rxn_id)

# dG0_data = data.getFreeEnergyData(GEM,
#                                   work_directory=par.work_directory,
#                                   pH_i=par.pH_i,
#                                   Ionic_strength=par.Ionic_strength,
#                                   dG0_uncertainty_threshold=par.alpha)

GEM
```

    Downloading package metadata...
    Fragments already downloaded
    Downloading package metadata...
    Fragments already downloaded
    Using license file C:\Users\tinta\gurobi.lic
    Academic license - for non-commercial use only
    Read LP format model from file C:\Users\tinta\AppData\Local\Temp\tmp9qizxd21.lp
    Reading time = 0.03 seconds
    : 1147 rows, 4276 columns, 17168 nonzeros
    





<table>
    <tr>
        <td><strong>Name</strong></td>
        <td>iML1515</td>
    </tr><tr>
        <td><strong>Memory address</strong></td>
        <td>0x0219cfa7f3c8</td>
    </tr><tr>
        <td><strong>Number of metabolites</strong></td>
        <td>1147</td>
    </tr><tr>
        <td><strong>Number of reactions</strong></td>
        <td>2138</td>
    </tr><tr>
        <td><strong>Number of groups</strong></td>
        <td>0</td>
    </tr><tr>
        <td><strong>Objective expression</strong></td>
        <td>1.0*BIOMASS_Ec_iML1515_core_75p37M - 1.0*BIOMASS_Ec_iML1515_core_75p37M_reverse_35685</td>
    </tr><tr>
        <td><strong>Compartments</strong></td>
        <td>cytosol, extracellular space, periplasm</td>
    </tr>
  </table>



