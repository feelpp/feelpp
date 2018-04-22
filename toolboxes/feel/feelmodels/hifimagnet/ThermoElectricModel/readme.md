# ThermoElectric Model

## Model

The models currently available are:

- thermoelectric-linear
- thermoelectric-nonlinear

## Boundary conditions

The fields are:



The condition types are:



## Options

All options are prefixed by `thermoelectric`

- order: approximation order (int)
- eps-coeff: penalisation parameter for regularized formulation (double)
- tolerance: tolerance for non linear computatin (ferromagnetism) (double)
- weakdir: use weakdir Dirichlet condition (bool)
- penaldir: penalisation parameter for the weak boudary condition (double)
- print_info: print some info (bool)
- export: export results (TODO: use model) (bool)
- model_json: json file containing the model (string)
