"""
Perspectives:
 . Other stencils using eigenvalues or t-designs
 . parallelization using MPI
 . riddance of used observations using pop method 
"""

import numpy as np
import scipy as sci
from scipy import linalg
import time

VERBOSE = not True
seconds = time.time()

""" This part contains basic functions not intrinsically related
    to the filter implementation 
"""

def tic(message = '', verbose: bool = False):
    """ Equivalent to tic in cpp """
    global seconds
    seconds = time.time()
    if VERBOSE or verbose:
        print(message)

def toc(message = '', verbose: bool = False):
    """ Equivalent to toc in cpp """
    if VERBOSE or verbose:
        print(message + ' {} seconds'.format(time.time() - seconds))

def isempty(liste: list) -> bool:
    return not bool(len(liste))

def isnumber(variable) -> bool:
    return type(variable) is float or type(variable) is int or type(variable) is np.float64

def balancedpartition(nb_data,nb_procs):
    """ For parallelization purpose using MPI
    
    Computes the distribution of data among the processors 
    as uniformly as possible
    """
    if nb_data % nb_procs == 0:
        return nb_procs*[nb_data//nb_procs]
    return (nb_procs-1)*[nb_data//(nb_procs-1)] + [nb_data%(nb_procs-1)]

def displacements(partition):
    """ For parallelization purpose using MPI

    Computes the index shift accordingly to the distribution in argument
    """
    tic()
    shift = [0]
    for i in range(1,len(partition)):
        shift.append(shift[i-1] + partition[i-1])
    toc()
    return shift


""" Here starts the implementation of the filter """

def magnitude_order(array):
    """ Returns the order of magnitude of a number or a vector.

    For the vector the L2 norm is used.
    """

    return 10**np.floor(np.log10(np.linalg.norm(array)))

def dispersion_matrix(array, optional_array = None):
    """ Computes the dispersion matrix associated with one or two vectors.

    input
        array:          np.array of size n
        optional_array: np.array of size m
    output
        return:         np.array of size (n,m)
    """

    if optional_array is None:
        return np.outer(array,array)
    else:
        return np.outer(array, optional_array)

def weighted_sum(element_list, weights = None):
    """ Computes the weighted sum of a list of elements.

    input
        element_list: list of objects for which + and * are defined
        weights:      list of floats of same length
    output
        return:       object of same type as any element
    """

    if weights is None:
        weights = np.ones(len(element_list))/len(element_list)
    elif type(weights) is list and len(weights) != len(element_list):
        raise ValueError('There must be the same number of states and weights')
    if type(element_list[0]) is State:
        result = State(dim = element_list[0].get_dim()) # sets zero State
        for state, weight in zip(element_list,weights):
            result += weight * state
    else:
        result = 0
        for element, weight in zip(element_list,weights):
            result += weight * element
    return result

class State:
    """ A State is the representation of a configuration 
    of the system of interest. The dimension can be accessed as an attribute.
    
    This class is used for States and for Observations since only their
    respective dimention may differ.

    input
        input: State or np.array or list of numbers or single number
        dim:   integer > 0
    
    attributes
        _values: State's coords in state space
        _dim:    state space's dimension

    methods
        __init__(input: list or np.array or number, dim: int) -> None
        set_values(values: list or np.array or number) -> None
        set_dim(integer: int) -> None
        get_values() -> np.ndarray
        get_dim() -> int
        round(order: int) -> None
        __add__(other: State) -> State
        __iadd__(other: State) -> State
        __sub__(other: State) -> State
        __mul__(other: State) -> State
        __rmul__(factor: float) -> State
        __truediv__(divider: float) -> State
        __matmul__(matrix) -> State
        __eq__(other) -> bool
        __getitem__(key) -> float
        __str__() -> None
    """
    def __init__(self, 
                 input = None, 
                 dim = 1):
        self._values = None
        self._dim = None
        if input is not None:
            if type(input) is State:
                self = input
            else:
                self.set_values(input)
        else:
            self.set_dim(dim)
        
    def set_values(self, values):
        """ Sets array as the state and computes the dimension 
        
        input
            values: np.array or list of numbers or single number
        """

        if type(values) is not np.ndarray and type(values) is not list and not isnumber(values):
            raise TypeError('state values must be 1D numpy.array or list or number')
        self._values = np.array([values])
        self._values.shape = (max(self._values.shape),)
        self._dim = len(self._values)

    def set_dim(self, integer: int):
        """ Sets the dimension of the state to integer 
        and sets the associated set to an array of zeros of that length 

        input
            integer > 0
        """
        if integer < 1:
            raise ValueError('the dimension must be at least 1')
        if type(integer) is not int:
            raise ValueError('the dimension must be an integer')
        self._dim = integer
        self._values = np.zeros((integer))

    def get_values(self) -> np.ndarray:
        """ Getter for values """
        return self._values

    def get_dim(self) -> int:
        """ Getter for dimension """
        return self._dim

    def round(self, order: int = 0) -> None:
        """ Method for rounding the State's values """

        return State(self.get_values().round(order))

    def __add__(self, other):
        if type(other) is State:
            return State(self.get_values() + other.get_values())
        elif isnumber(other):
            return State(self.get_values() + other)
        else:
            raise TypeError('invalid sum term')

    def __iadd__(self, other):
        if type(other) is State:
            self._values += other.get_values()
        elif isnumber(other):
            self._values += other
        else:
            raise TypeError('invalid increment')
        self._values.shape = (max(self._values.shape),)
        return self

    def __mul__(self, factor: float):
        return State(factor * self.get_values())

    def __sub__(self, other):
        if type(other) is State:
            return self + (-1)*other
        elif isnumber(other):
            return State(self.get_values() - other)
        else:
            raise TypeError('invalid substraction term')

    def __rmul__(self, factor: float):
        return self * factor

    def __truediv__(self, divider: float):
        return self*(1/divider)

    def __matmul__(self, matrix):
        return State(np.array(matrix @ self.get_values()))

    def __eq__(self, other):
        return np.array_equal(self.get_values(), other.get_values()) and self.get_dim() == other.get_dim()

    def __getitem__(self, key):
        return self.get_values()[key]

    def __str__(self):
        return "State of dimension {}".format(self.get_dim())


class EnsembleTools:
    """ Contains all the methods and parameters used in computations related
    to the ensemble

    ensemble: iterable list of states
    sigma_scheme is an iterable list of shift states (default is UKF friendly)

    input
        factor:      float, factor used in the computation of sigma-points
        main_weight: float, weight of the central estimation
        dim:         dictionary of int with entries state and measure

    attributes
        ensemble:    list of States
        size:        depends on the dimension and on the method (default is UKF)
        dim:         dictionary of int with entries state and measure
        factor:      float, factor used in the computation of sigma-points
        covariances: dictionary of arrays with entries state, observation, cross and measure

    methods
        __init__(factor: float, main_weight: float, dim: dict) -> None
        set_dim(value: int, type: str) -> None
        set_weights(value: float) -> None
        set_covarianves(value: np.array or float, type: str) -> None
        produce_ensemble(state: State) -> None
    """

    def __init__(self, 
                 factor: float, 
                 main_weight: float, 
                 dim: dict = {'state': 1, 'measure': 1}):
        self.ensemble = None
        self.size = 2*dim['state']+1
        self.dim = dim
        self.factor = factor
        self.set_weights(main_weight)
        self.covariances = {}
        for type in ['state', 'observation', 'cross', 'measure']:
            self.set_covariance(1, type)

    def set_dim(self, value: int, type: str):
        if type in ['state','measure']:
            self.dim[type] = value
            self.size = 2*self.dim['state']+1
            self.set_weights(self.weights[0])
            self.set_covariance(1, type)
        else:
            raise ArgumentError('\'type\' must be \'state\' or \'measure\'')

    def set_weights(self, value: float):
        """ Computes the weights associated with the stencil size 
        and given central weight.
        """

        self.weights = [value] + (self.size-1)*[(1-value)/(2*self.dim['state'])]

    def set_covariance(self, value, type: str):
        """ Sets the state or observation covariance. 

        If a number is provided, the resulting matrix is number*Id.
        Else a square matrix of appropriate size is required.
        Zero is not allowed since it would mean maximal confidence.
        """
        if type not in ['state', 'observation', 'cross', 'measure']:
            raise ValueError('Covariance type is either \'state\', \'observation\', \'cross\' or \'measure\'')
        if np.linalg.norm(value) != 0:
            if isnumber(value):
                if type is 'measure' or type is 'state':
                    self.covariances[type] = value * np.eye(self.dim[type])
            else:
                self.covariances[type] = value
        else:
            raise ValueError('A covariance cannot be zero')

    def produce_ensemble(self, state: State):
        """ Computes the stencil of sigma-points 
        using the matrix square root method.

        Only UKF is implemented.
        """

        self.ensemble = [state]
#        shift_matrix = np.linalg.eig(0.1*self.covariances['state'])[1]
        shift_matrix = sci.linalg.sqrtm(self.factor*self.covariances['state'])
#        print('\n shift2:\n', (self.factor*self.covariances['state']).round(4).reshape(1,-1).tolist()[0])
        for column in range(len(shift_matrix)):
            self.ensemble.append(state + State(shift_matrix[:,column]))
            self.ensemble.append(state - State(shift_matrix[:,column]))

class Filter:
    """ A Filter object is the data of:
    . a list of estimated State
    . a list of observed State
    . some parameters

    The ensemble size is:
    . 1 for simple Kalman filter (not implemented)
    . any for ensemble Kalman filter (not implemented)
    . 2*dim+1 for unscented Kalman filter (default)

    input
        dimension:      dictionary with state and measure entries, both >0 integers
        main_weight:    weight associated with the central estimation
        initial_state:  starting point
        forecast_state: function mapping a state onto the predicted state one time step further
        forecast_obs:   function mapping a state onto associated expected measurement

    attributes
        forecast_state:    function mapping a state onto the predicted state one time step further
        forecast_obs:      function mapping a state onto associated expected measurement
        real_observations: list of provided measurements
        gain:              current Kalman gain
        ts:                current timestep
        max_ts:            last timestep
        tools:             EnsembleTools instance
        analyzed_states:   time series of best-knowledge estimated states

    methods
        __init__(dimension: dict, main_weight: float, forecast_state: func, forecast_obs: func)
        set_state(state: State) -> None
        load_measurements(data: list) -> None
        set_forecast_function(function: func) -> None
        set_observation_function(function: func) -> None
        get_last_state() -> State
        get_last_obs() -> State
        get_ts() -> int
        step_ts() -> None
        compute_gain(x_f, x_dag, y_dag) -> None
        analyze(x_f, x_dag, y_dag) -> None
        forecast() -> None
        filter(initial_guess: State, data: list) -> None
        extract_analyzed_states(components: str or list of int) -> np.array
        __str__() -> None
    """

    def __init__(self, 
                 dimension: dict = {'state': 1, 'measure': 1}, 
                 main_weight: float = 0.5,
                 initial_state: State = None,
                 forecast_state = lambda x : x, 
                 forecast_obs = lambda x : x):

        self.forecast_state = forecast_state
        self.forecast_obs = forecast_obs
        
        self.real_observations = []
        self.gain = None
        self.ts = 0
        self.max_ts = None
        self.tools = EnsembleTools(factor = dimension['state']/(1-main_weight), 
                                   main_weight = main_weight,
                                   dim = dimension)
        self.analyzed_states = []
        if initial_state is not None:
            self.set_state(initial_state)
        else:
            self.set_state(State(dim = dimension['state']))
        self.estimate_observations = [self.forecast_obs(self.get_last_state())]

    def set_state(self, state: State):
        if type(state) is State:
            self.analyzed_states.append(state)
        else:
            self.analyzed_states.append(State(state))
        self.tools = EnsembleTools(factor = self.get_last_state().get_dim()/(1-self.tools.weights[0]),
                                   main_weight = self.tools.weights[0])
        self.tools.set_dim(self.get_last_state().get_dim(), 'state')

    def load_measurements(self, list):
        tic()
        self.real_observations = []
        for state in list:
            if type(state) is State:
                self.real_observations.append(state)
            else:
                self.real_observations.append(State(input = state))
        self.max_ts = len(self.real_observations)
        self.tools.set_dim(self.real_observations[0].get_dim(), 'measure')
        toc('{} observations loaded in'.format(self.max_ts))

    def set_forecast_function(self, function):
        self.forecast_state = function
    
    def set_observation_function(self, function):
        self.forecast_obs = function
        
    def get_last_state(self) -> State:
        return self.analyzed_states[-1]

    def get_last_obs(self) -> State:
        return self.estimate_observations[-1]

    def get_ts(self) -> int:
        return self.ts

    def step_ts(self):
        self.ts += 1

    def compute_gain(self, x_f, x_dag, y_dag):
        """ Computes the gain """

        self.tools.set_covariance(
            weighted_sum(
                element_list = [dispersion_matrix(x.get_values()-x_f.get_values(), y.get_values() - self.get_last_obs().get_values()) for x,y in zip(x_dag, y_dag)],
                weights = self.tools.weights
                ),
            'cross'
            )
        self.tools.set_covariance(
            weighted_sum(
                element_list = [dispersion_matrix(y.get_values() - self.get_last_obs().get_values()) for y in y_dag],
                weights = self.tools.weights
                ),
            'observation'
            )
        self.gain = self.tools.covariances['cross'] \
                    @ np.linalg.inv(self.tools.covariances['observation'] \
                                    + self.tools.covariances['measure'])

    def analyze(self, x_f, x_dag, y_dag):
        """ Computes the analysis formula using most recent data 
        and provided forecast state
        """

        self.compute_gain(x_f, x_dag, y_dag)
        self.analyzed_states.append(x_f + ((self.real_observations[self.ts] - self.get_last_obs()) \
                                            @ self.gain))

    def forecast(self):
        """ The first part computes the sigma-points and transforms them.
        Then the associated preficted observations are computed.

        The second part uses these informations along with a real observation
        to compute the analyzed (i.e. best knowledge) state
        """

        tic()
        if isempty(self.real_observations):
            raise ValueError('Obervation data missing')
        if self.ts >= self.max_ts:
            raise RecursionError('Last time step exceeded')
        
        ensemble_state = []
        ensemble_obs = []

        self.tools.produce_ensemble(self.get_last_state())
        for guess in self.tools.ensemble:
            sigma_point = self.forecast_state(guess)
            ensemble_state.append(sigma_point)
            ensemble_obs.append(self.forecast_obs(sigma_point))
        self.tools.ensemble = ensemble_state

        forecast_state = weighted_sum(
            element_list = ensemble_state,
            weights = self.tools.weights
            )
        self.estimate_observations.append(
            weighted_sum(
                element_list = ensemble_obs,
                weights = self.tools.weights
                )
            )

        self.analyze(forecast_state, ensemble_state, ensemble_obs)
        self.tools.set_covariance(
            (np.eye(self.tools.dim['state']) - self.gain @ self.tools.covariances['measure']) @ \
            weighted_sum(
                element_list = [dispersion_matrix(s.get_values() - self.get_last_state().get_values()) for s in self.tools.ensemble],
                weights = self.tools.weights
                ),
                type = 'state'
        )    
        self.step_ts()

        del sigma_point, ensemble_state, ensemble_obs, forecast_state
        toc('Step {} achieved in'.format(self.get_ts()))

    def filter(self, initial_guess = None, data = None):
        if isempty(self.analyzed_states) and initial_guess is None:
            raise ValueError('Please provide an initial state')
        if isempty(self.real_observations) and data is None:
            raise ValueError('Please provide data')
        if initial_guess is not None:
            self.__init__(initial_state = initial_guess)
        if data is not None:
            self.load_measurements(data)

        while self.ts < self.max_ts:
            self.forecast()

    def extract_analyzed_states(self, components = 'all'):
        """ Extracts analyzed states trajectories.

        The argument \"component\" is either \"all\" or a list of indices.
        """

        output = []
        indices = range(self.tools.dim['state']) if components == "all" else components
        for state in self.analyzed_states:
            output.append(state.get_values()[indices].reshape((len(indices),)))
        return np.array(output[1:]).squeeze()

    def __str__(self):
        message = "Filter at time step {}\n".format(self.get_ts()) \
            + "    State dimension = {}\n".format(self.tools.dim['state'])
        if isempty(self.real_observations):
            message += "    no observations"
        else:
            message += "    Obs   dimension = {}".format(self.tools.dim['measure'])
        return message
