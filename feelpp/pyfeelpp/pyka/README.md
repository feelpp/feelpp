# PyKa module

## Public classes

### State

A State is the representation of a configuration of the system of interest. The dimension can be accessed as an attribute. This class is used for both States and Observations since only their respective dimention may differ.

* input
    * input: State or np.array or list of numbers or single number
    * dim:   integer > 0
    
* attributes
    * _values: State's coords in state space
    * _dim:    state space's dimension

* methods
    * __init__(input: list or np.array or number, dim: int) -> None
    * set_values(values: list or np.array or number) -> None
    * set_dim(integer: int) -> None
    * get_values() -> np.ndarray
    * get_dim() -> int
    * round(order: int) -> None
    * __add__(other: State) -> State
    * __iadd__(other: State) -> State
    * __sub__(other: State) -> State
    * __mul__(other: State) -> State
    * __rmul__(factor: float) -> State
    * __truediv__(divider: float) -> State
    * __matmul__(matrix) -> State
    * __eq__(other) -> bool
    * __getitem__(key) -> float
    * __str__() -> None

### Filter

A Filter object is the data of:
* a list of estimated State
* a list of observed State
* some parameters

The ensemble size is:
* 1 for simple Kalman filter (not implemented)
* any for ensemble Kalman filter (not implemented)
* 2*dim+1 for unscented Kalman filter (default)


* input
    * dimension:      dictionary with state and measure entries, both >0 integers
    * main_weight:    weight associated with the central estimation
    * initial_state:  starting point
    * forecast_state: function mapping a state onto the predicted state one time step further
    * forecast_obs:   function mapping a state onto associated expected measurement

* attributes
    * forecast_state:    function mapping a state onto the predicted state one time step further
    * forecast_obs:      function mapping a state onto associated expected measurement
    * real_observations: list of provided measurements
    * gain:              current Kalman gain
    * ts:                current timestep
    * max_ts:            last timestep
    * tools:             EnsembleTools instance
    * analyzed_states:   time series of best-knowledge estimated states

* methods
    * __init__(dimension: dict, main_weight: float, forecast_state: func, forecast_obs: func)
    * set_state(state: State) -> None
    * load_measurements(data: list) -> None
    * set_forecast_function(function: func) -> None
    * set_observation_function(function: func) -> None
    * get_last_state() -> State
    * get_last_obs() -> State
    * get_ts() -> int
    * step_ts() -> None
    * compute_gain(x_f, x_dag, y_dag) -> None
    * analyze(x_f, x_dag, y_dag) -> None
    * forecast() -> None
    * filter(initial_guess: State, data: list) -> None
    * extract_analyzed_states(components: str or list of int) -> np.array
    * __str__() -> None

## Private class

### EnsembleTools

Contains all the methods and parameters used in computations related to the ensemble.

* input
    * factor:      float, factor used in the computation of sigma-points
    * main_weight: float, weight of the central estimation
    * dim:         dictionary of int with entries state and measure

* attributes
    * ensemble:    list of States
    * size:        depends on the dimension and on the method (default is UKF)
    * dim:         dictionary of int with entries state and measure
    * factor:      float, factor used in the computation of sigma-points
    * covariances: dictionary of arrays with entries state, observation, cross and measure

* methods
    * __init__(factor: float, main_weight: float, dim: dict) -> None
    * set_dim(value: int, type: str) -> None
    * set_weights(value: float) -> None
    * set_covarianves(value: np.array or float, type: str) -> None
    * produce_ensemble(state: State) -> None

## Pytest

### State class tests

* test_init: state initialization
* test_addition: state addition
* test_multiplication: scalar multiplication
* test_weighted_sum: state weighted sum

### Filter class tests

* test_init: filter initialization
* test_load_data: data loading
* test_simple_direct_filter: direct 1D filter with identity forecasts
* test_multiple_direct_filter: direct nD filter with identity forecasts
* test_sum_diff: indirect measurement
    * real hidden state is [1,1]
    * observation is [sum,difference] with small uncertainty
* test_noisy_tracking: noisy sine. Noisy data is filtered using sine's order 3 Taylor expansion.
* test_constant_velocity_estimation: constant velocity particle
    * Noisy position is measured
    * Trajectory is filtered and velocity is estimated