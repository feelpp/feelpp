"""
Perspectives:
 . parallelization using MPI
 . riddance of used observations using pop method 
"""

import numpy as np
import scipy as sci
import scipy.linalg
import time

VERBOSE = True
seconds = time.time()

""" This part contains basic functions not intrinsically related
    to the filter implementation 
"""

def tic(message = ''):
    """ Equivalent to tic in cpp """
    global seconds
    seconds = time.time()
    if VERBOSE:
        print(message)

def toc(message = ''):
    """ Equivalent to toc in cpp """
    if not VERBOSE and message is not '':
        print(message + ' {} seconds'.format(time.time() - seconds, '1.2'))
    if VERBOSE:
        print(message + ' {} seconds'.format(time.time() - seconds))

def isempty(liste: list) -> bool:
    return not bool(len(liste))

def isnumber(variable) -> bool:
    return type(variable) is float or type(variable) is int

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

def weighted_sum(element_list, weights = None):
    if weights is None:
        weights = np.ones(len(element_list))/len(element_list)
    elif type(weights) is list and len(weights) != len(element_list):
        raise ValueError('There must be the same number of states and weights')
    if type(element_list[0]) is State:
        result = State(dim = element_list[0].get_dim())
        for state, weight in zip(element_list,weights):
            result += weight * state
    else:
        result = 0
        for element, weight in zip(element_list,weights):
            result += weight * element
    return result

class State:
    """ A State is the numpy.array representation of a configuration 
    of the system of interest. The dimension can be accessed as an attribute.
    
    This class is used for States and for Observations since only their
    respective dimention may differ.
    """
    def __init__(self, input = None, dim = 1):
        self._values = None
        self._dim = None
        if input is not None:
            self.set_values(input)
        else:
            self.set_dim(dim)
        
    def set_values(self, values):
        """ Sets array as the state and computes the dimension """
        if type(values) is not np.ndarray and type(values) is not list and not isnumber(values):
            raise TypeError('state description must be 1D numpy.array or list or number')
        self._values = np.array([values])
        self._values.shape = (max(self._values.shape),)
        self._dim = len(self._values)

    def set_dim(self, integer: int):
        """ Sets the dimension of the state to integer 
        and sets the associated set to an array of zeros of that length 
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

    def __str__(self):
        return "State of dimension {}".format(self.get_dim())


class EnsembleTools:
    """ Contains all the methods and parameters used in computations related
    to the ensemble

    ensemble is an iterable list of states
    sigma_scheme is an iterable list of shift states (default is UKF friendly)
    """

    def __init__(self, 
                 factor: float, 
                 main_weight: float, 
                 dim: int):
        self.ensemble = None
        self.size = 2*dim+1
        self.dim = dim
        self.factor = factor
        self.weights = self.set_weights(main_weight, dim)
        self.covariances = {'state': None, 'observation': None}
        for type in ['state', 'observation', 'cross', 'measure']:
            self.set_covariance(1, type)

    def set_weights(self, value: float, dim: int):
        """ Computes the weights associated with the stencil size 
        and given central weight.
        """

        return [value] + (self.size-1)*[(1-value)/(2*self.dim)]

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
                self.covariances[type] = value * np.eye(self.dim)
            else:
                self.covariances[type] = value

    def produce_ensemble(self, state: State):
        """ Computes the stencil of sigma-points 
        using the matrix square root method.

        Only UKF fashion is implemented.
        """

        self.ensemble = [state]
        shift_matrix = np.sqrt(self.factor) * sci.linalg.sqrtm(self.covariances['state'])
        print(self.covariances['state'])
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

    forecast_state and forecast_obs are functions mapping a state 
    to the predicted state one time step further ; default is identity
    """

    def __init__(self, 
                 dimension: int = 1, 
                 main_weight: float = 0.5,
                 initial_state: State = None,
                 forecast_state = lambda x : x, 
                 forecast_obs = lambda x : x):

        self.forecast_state = forecast_state
        self.forecast_obs = forecast_obs
        self.estimate_states = [State(input = initial_state)] if initial_state is not None else [State(dim = dimension)]
        self.estimate_observations = [self.forecast_obs(self.get_last_state())]
        self.real_observations = None
        self.gain = None
        self.ts = 0
        self.max_ts = None
        self.tools = EnsembleTools(factor = dimension/(1-main_weight), 
                                   main_weight = main_weight,
                                   dim = dimension)

    def set_state(self, state: State):
        if type(state) is State:
            self.estimate_states.append(state)
        else:
            self.estimate_states.append(State(state))

    def load_real_observations(self, list):
        tic()
        self.real_observations = []
        for state in list:
            if type(state) is State:
                self.real_observations.append(state)
            else:
                self.real_observations.append(State(input = state))
        self.max_ts = len(self.real_observations)
        toc('Observations loaded in')

    def set_forecast_function(self, function):
        self.forecast_state = function
    
    def set_observation_function(self, function):
        self.forecast_obs = function
        
    def get_last_state(self) -> State:
        return self.estimate_states[-1]

    def get_last_obs(self) -> State:
        return self.estimate_observations[-1]

    def get_ts(self) -> int:
        return self.ts

    def step_ts(self):
        self.ts += 1
        if self.max_ts is not None:
            if self.ts == self.max_ts:
                print('Last time step is reached.')

    def compute_gain(self, x_f, x_dag, y_dag):
        """ Computes the gain (transposed for easier further computation) """
        S_xy = []
        S_yy = []
        for x, y in zip(x_dag, y_dag):
            S_xy.append(np.mat(x.get_values() - x_f.get_values()).T \
                               @ np.mat(y.get_values() - self.get_last_obs().get_values()))
            S_yy.append(np.mat(y.get_values() - self.get_last_obs().get_values()).T \
                        @ np.mat(y.get_values() - self.get_last_obs().get_values()))
        self.tools.set_covariance(weighted_sum(S_xy, weights = self.tools.weights), 'cross')
        self.tools.set_covariance(weighted_sum(S_yy, weights = self.tools.weights), 'observation')
        self.gain = self.tools.covariances['cross'] \
                    @ np.linalg.inv(self.tools.covariances['observation'] \
                                    + self.tools.covariances['measure'])

        del S_xy, S_yy

    def analyze(self, x_f, x_dag, y_dag):
        """ Computes the analysis formula using most recent data 
        and provided forecast state
        """

        self.compute_gain(x_f, x_dag, y_dag)
        self.estimate_states.append(x_f + ((self.real_observations[self.ts] - self.get_last_obs()) \
                                            @ self.gain))

    def forecast(self):
        """ The first part computes the sigma-points and transforms them.
        Then the associated preficted observations are computed.

        The second part uses these informations along with a real observation
        to compute the analyzed (i.e. best knowledge) state
        """

        tic()
        if self.real_observations is None:
            raise ValueError('Obervation data missing')
        if self.ts == self.max_ts:
            raise RecursionError('Last time step reached')
        
        ensemble_state = []
        ensemble_obs = []

        self.tools.set_covariance(np.outer(self.get_last_state().get_values(),
                                           self.get_last_state().get_values()),
                                  type = 'state')
        self.tools.set_covariance(np.outer(self.real_observations[self.ts].get_values(),
                                           self.real_observations[self.ts].get_values()),
                                  type = 'observation')

        self.tools.produce_ensemble(self.get_last_state())
        for guess in self.tools.ensemble:
            sigma_point = self.forecast_state(guess)
            ensemble_state.append(sigma_point)
            ensemble_obs.append(self.forecast_obs(sigma_point))
            #print(sigma_point.get_values())
            print(guess.get_values())

        forecast_state = weighted_sum(element_list = ensemble_state,
                                      weights = self.tools.weights)
        self.estimate_observations.append(weighted_sum(element_list = ensemble_obs,
                                                       weights = self.tools.weights))

        self.analyze(forecast_state, ensemble_state, ensemble_obs)
        self.step_ts()

        del sigma_point, ensemble_state, ensemble_obs, forecast_state
        toc('Step {} achieved in'.format(self.get_ts()))

    def filter(self, initial_guess, data):
        pass

    def __str__(self):
        message = "Filter at time step {}\n".format(self.get_ts()) \
            + "    State dimension = {}\n".format(self.get_last_state().get_dim())
        if self.real_observations is None:
            message += "    no observations"
        else:
            message += "    Obs   dimension = {}".format(self.get_last_obs().get_dim())
        return message
