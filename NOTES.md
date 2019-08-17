
## Fixing certain concentrations
In Thermodynamics-based metabolic analysis, authors fixed the concentration of intracellular protons to 1e-7M, which corresponds to pH = 7. Also, internal oxygen concentration cannot be higher than that of outside, since no reaction in E. coli produces oxygen and oxygen diffuses passively following a gradient. Update: already implemented.

## Using a pH range instead of pH = 7
Perhaps good idea to use a range of pH values, such as 7.2 - 7.7 when calculating $\Delta G$s, one would use then the minimum and maximum values of these in the bounds. Update: already implemented.

## UPDATE: 02/VI/2019
Tried the code with 1% deviation of the dG0 mean to filter the data (Cond1). Obtained a graph with only 5 metabolites and took around 5000 seconds. Clearly need to do preprocessing if I want to soften the filter. I could apply the same idea used in the flux orders: sample the space of possible concentrations and discard the pairs with conflicting order relation.

Should also organize the code: it is a mess!

## Organizing the code
After reading and thinking about it, I believe now that the best way to organize code in python is through modules, as opposed to classes. The only tricky thing with python is that namespaces are not shared between modules. So for instance, if I have a module named parameters and I want to use some of the variables defined there in another module called module1, then I would have to do something like this in module1

```python3
import parameters as p
p.parameter1
```

or alternatively

```python3
from parameters import parameter1, parameter2
parameter1 + parameter2
```

These two alternatives would work to transfer objects between modules.

## Reducing computational time
Sampling takes too much time when delta_dG0 > 0.1 because there are too many constraints. What ca we do?

1. Increase tolerance of the objective function: I do no care if the optimum is reached since I just want to generate a random point inside the feasible region.
2. Use the MILP to obtain the path of alternative suboptimal solutions as my sample. Perhaps glpk even provides this option.
3. Use Gurobi
4. Perhaps there is another, cheaper, method to sample the feasible space via a MILP? i.e., a different objective function that is less demanding.

Update: I implemented some of these ideas and I have been able to reduce the computational time significantly.

## A coment on results so far
I've seen that the number of ordered pairs is extremely sensitivy to:
1. The global concentration lower bound (almost no pairs when x_min < 1e-4)
2. The minimum pH value (many pairs when pH_min > 7.1)

Why is that? What can we conclude from this? What other experiments, observations are relevant?
