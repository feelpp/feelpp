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
    if VERBOSE:
        print(message + str(time.time() - seconds) + ' seconds')

def isempty(liste: list) -> bool:
    return not bool(len(liste))

def isnumber(variable) -> bool:
    return type(variable) is float or type(variable) is int

def balancedpartition(nb_data,nb_procs):
    """ For parallelization purpose using MPI
    
    Computes the distribution of data among the processors 
    as evenly as possible
    """
    tic()
    partition = [] #np.zeros(nb_procs, dtype=np.int8)
    for i in range(nb_procs):
        nb_partitions = int(round(nb_data/(nb_procs-i)))
        partition.append(nb_partitions)
        nb_data -= nb_partitions
    toc()
    return partition

def displacements(partition):
    """ For parallelization purpose using MPI

    Computes the index shift accordingly to the distribution in argument
    """
    tic()
    shift = []
    for i in range(1,len(partition)):
        shift.append(displacement[i-1] + partition[i-1])
    toc()
    return shift


""" Here starts the implementation of the filter """

def magnitude_order(array):
    return 10**np.floor(np.log10(np.linalg.norm(array)))

def weighted_sum(element_list, weights = None):
    if weights is None:
        weights = np.ones(len(element_list))/len(element_list)
    elif type(weights) is list and len(weights) != len(element_list):
        raise ValueError('There must be the same number of states and weights')
    if type(element_list[0]) is State:
        result = State(element_list[0].get_dim())
        for state, weight in zip(element_list,weights):
            result += state * weight
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
    def __init__(self, input):
        self._values = None
        self._dim = None
        if type(input) is int:
            self.set_dim(input)
        else:
            self.set_values(input)

    def set_values(self, values):
        """ Sets array as the state and computes the dimension """
        if type(values) is not np.ndarray and type(values) is not list:
            raise TypeError('state must be 1D np.array or list')
        self._values = np.array(values)
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
        self.factor = factor
        self.weights = self.set_weights(main_weight, dim)
        self.state_cov = self.set_state_cov(1)
        self.obs_cov = self.set_state_cov(1)

    def set_weights(self, value: float, dim: int):
        return [value] + (self.size-1)*[(1-value)/(2*dim)]

    def set_state_cov(self, value):
        if isnumber(value):
            self.state_cov = value * np.eye(dim)
        else:
            self.state_cov = value

    def set_obs_cov(self, value):
        if isnumber(value):
            self.obs_cov = value * np.eye(dim)
        else:
            self.obs_cov = value

    def produce_ensemble(self, state: State):
        self.ensemble = [state]
        shift_matrix = np.sqrt(self.factor) * sci.linalg.sqrtm(self.state_cov)
        signs = [0] + state.get_dim()*[1] + state.get_dim()*[-1]
        for column in range(len(shift_matrix)):
            self.ensemble.append(state + State(shift_matrix[:,column]))
            self.ensemble.append(state - State(shift_matrix[:,column]))

class Filter:
    """ A Filter object is the data of:
    . estimated State
    . observed State

    The ensemble size is:
    . 1 for simple Kalman filter
    . any for ensemble Kalman filter
    . 2*dim+1 for unscented Kalman filter (default)

    forecast_state and forecast_obs are functions mapping a state 
    to the predicted state one time step further ; default is identity
    """

    def __init__(self, 
                 dim: int = 1, 
                 main_weight: float = 0.5,
                 forecast_state = lambda x : x, 
                 forecast_obs = lambda x : x):

        self.estimate_states = [State(dim)]
        self.estimate_observations = []
        self.real_observations = None
        self.gain = None
        self.ts = 0
        self.tools = EnsembleTools(factor = dim/(1-main_weight), 
                                   main_weight = main_weight,
                                   dim = dim)
        self.forecast_state = forecast_state
        self.forecast_obs = forecast_obs

    def set_state(self, state: State):
        if type(state) is State:
            self.estimate_states.append(state)
        else:
            self.estimate_states.append(State(state))

    def load_real_observations(self, list):
        self.real_observations = []
        for state in list:
            if type(state) is State:
                self.real_observations.append(state)
            else:
                self.real_observations.append(State(state))

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

    def compute_gain(self, x_f, x_dag, y_dag):
        """ Computes the transposed gain """
        S_xy = []
        S_yy = []
        for x, y in zip(x_dag, y_dag):
            S_xy.append(np.mat(x.get_values() - x_f.get_values()).T \
                               @ np.mat(y.get_values() - self.get_last_obs().get_values()))
            S_yy.append(np.mat(y.get_values() - self.get_last_obs().get_values()).T \
                        @ np.mat(y.get_values() - self.get_last_obs().get_values()))
    
        self.gain = weighted_sum(S_xy, weights = self.tools.weights) \
               @ np.linalg.inv(weighted_sum(S_yy, weights = self.tools.weights) \
                               + sci.sqrtm(self.tools.obs_cov))

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

        if self.real_observations is None:
            raise ValueError('Obervation data missing')
        
        ensemble_state = []
        ensemble_obs = []

        self.tools.set_state_cov(np.outer(self.get_last_state().get_values(),
                                          self.get_last_state().get_values()))
        self.tools.set_obs_cov(np.outer(self.real_observations[self.ts],
                                        self.real_observations[self.ts]))

        self.tools.produce_ensemble(self.get_last_state())
        for guess in self.tools.ensemble:
            sigma_point = self.forecast_state(guess)
            ensemble_state.append(sigma_point)
            ensemble_obs.append(self.forecast_obs(sigma_point))

        forecast_state = weighted_sum(ensemble_state,
                                      weights = self.tools.weights)
        self.estimate_observations.append(weighted_sum(ensemble_obs,
                                          weights = self.tools.weights))

        self.analyze(forecast_state, ensemble_state, ensemble_obs)
        self.step_ts()

    def filter(self, initial_guess, data):
        pass

    def __str__(self):
        message = "Filter at time step {}\n".format(self.get_ts()) \
            + "    State dimension = {}\n".format(self.get_last_state().get_dim())
        if self.observations is None:
            message += "    no observations"
        else:
            message += "    Obs   dimension = {}".format(self.get_last_obs().get_dim())
        return message
