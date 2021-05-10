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

def weighted_mean(list_of_states, weights = None):
    result = State(list_of_states[0].get_dim())
    if weights is None:
        weights = np.ones(len(list_of_states))/len(list_of_states)
    elif type(weights) is list and len(weights) != len(list_of_states):
        raise TypeError('There must be the same number of states and weights')
    for state, weight in zip(list_of_states,weights):
        result += state * weight
    return result

class State:
    """ A State is the numpy.array representation of a configuration 
    of the system of interest. The dimension can be accessed as an attribute.
    
    This class is used for States and for Observations since only their
    respective dimention may differ.
    """
    def __init__(self, input): #dimension: int = 1, values: list = [0]):
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
            raise TypeError('the dimension must be at least 1')
        if type(integer) is not int:
            raise TypeError('the dimension must be an integer')
        self._dim = integer
        self._values = np.zeros((integer))

    def get_values(self) -> np.ndarray:
        """ Getter for values """
        return self._values

    def get_dim(self) -> int:
        """ Getter for dimension """
        return self._dim

    def __add__(self, other):
        return State(self.get_values() + other.get_values())

    def __iadd__(self, other):
        self._values += other.get_values()
        return self

    def __mul__(self, factor: float):
        return State(factor * self.get_values())

    def __sub__(self, other):
        return self + (-1)*other

    def __rmul__(self, factor: float):
        return self * factor

    def __truediv__(self, divider: float):
        return self*(1/divider)

    def __str__(self):
        return "State of dimension {}".format(self.get_dim())


class EnsembleTools:
    """ Contains all the methods and parameters used in computations related
    to the ensemble

    ensemble is an iterable list of states
    sigma_scheme is an iterable list of shift states (default is UKF friendly)
    """

    def __init__(self, factor: float, main_weight: float, dim: int):
        self.ensemble = None
        self.size = 2*dim+1
        self.factor = factor
        self.weights = self.set_weights(main_weight, dim)
        self.state_cov = np.eye(dim)

    def set_weights(self, value: float, dim: int):
        return [value] + (self.size-1)*[(1-value)/(2*dim)]

    def produce_ensemble(self, state: State):
        self.ensemble = [state]
        shift_matrix = np.sqrt(self.factor) * sci.linalg.sqrtm(self.state_cov)
        signs = [0] + state.get_dim()*[1] + state.get_dim()*[-1]
        for column in range(len(shift_matrix)):
            self.ensemble.append(state + State(shift_matrix[:,column]))
            self.ensemble.append(state - State(shift_matrix[:,column]))
        return self.ensemble

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

        self.states = [State(dim)]
        self.observations = None
        self.ensemble = None
        self.ts = 0
        self.tools = EnsembleTools(factor = dim/(1-main_weight), 
                                   main_weight = main_weight,
                                   dim = dim)
        self.forecast_state = forecast_state
        self.forecast_obs = forecast_obs

    def set_state(self, state: State):
        self.states.append(state)

    def set_obs(self, obs: State):
        self.observations.append(obs)
        
    def get_last_state(self) -> State:
        return self.states[-1]

    def get_last_obs(self) -> State:
        return self.observations[-1]

    def get_ts(self) -> int:
        return self.ts

    def forecast(state: State):
        self.state = weighted_mean(self.tools.produce_ensemble(state),
                                   weights = self.tools.weights)

    def __str__(self):
        message = "Filter at time step {}\n".format(self.get_ts()) \
            + "    State dimension = {}\n".format(self.get_last_state().get_dim())
        if self.observations is None:
            message += "    no observations"
        else:
            message += "    Obs   dimension = {}".format(self.get_last_obs().get_dim())
        return message
