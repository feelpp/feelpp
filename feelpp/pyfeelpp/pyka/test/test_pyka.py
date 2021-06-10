from pyka import *
import numpy as np
import pytest

class TestState:

    def test_init(self):
        """ Test: state initialization """
        s = State(input = [1,2,3,4,5])
        assert s.get_dim() == 5
        assert np.array_equal(s.get_values(), np.array([1,2,3,4,5]))
        assert str(s) == "State of dimension 5"

    def test_addition(self):
        """ Test: state addition """
        assert State(input = [1,2]) + State(input = [2,3]) == State(input = [3,5])

    def test_multiplication(self):
        """ Test: scalar multiplication """
        assert 3 * State(input = [1,2,3]) == State(input = [3,6,9])

    def test_weighted_sum(self):
        """ Test: state weighted sum """
        list = [State(input = [1,2]), State(input = [2,3]), State(input = [3,-4])]
        weights = [2, -1, 1]
        assert weighted_sum(list, weights) == State(input = [3, -3])

class TestFilter():

    def test_init(self):
        """ Test: filter initialization """
        f = Filter(initial_state = [1,1,1])
        assert f.get_last_state().get_dim() == 3
        assert np.array_equal(f.get_last_state().get_values(), np.array([1,1,1]))
        assert str(f) == "Filter at time step 0\n    State dimension = 3\n    no observations"

    def test_load_data(self):
        """ Test: data loading """
        f = Filter()
        data = [np.array([1,2,3]),np.array([2,3,4]),np.array([3,4,5]),np.array([4,5,6]),np.array([5,6,7])]
        f.load_measurements(data)
        assert f.real_observations[0] == State(input = [1,2,3])
        assert f.max_ts == 5

    def test_simple_direct_filter(self):
        """ Test: direct 1D filter with identity forecasts """
        start = np.random.random()
        goal = np.random.random()
        f = Filter(initial_state = start)
        f.load_measurements(100*[goal])
        f.filter()
        assert np.linalg.norm(f.get_last_state().get_values() - State(input = goal).get_values()) < 0.05

    def test_multiple_direct_filter(self):
        """ Test: direct nD filter with identity forecasts """
        dim = 10
        start = np.random.random(dim)
        goal = np.random.random(dim)
        f = Filter(initial_state = start)
        f.load_measurements(100*[goal])
        f.filter()
        assert np.linalg.norm(f.get_last_state().get_values() - State(input = goal).get_values()) < 0.05

    def test_sum_diff(self):
        """ Test: indirect measurement
        
        - real hidden state is [1,1]
        - observation is [sum,difference] with small uncertainty
        """

        forecast_obs = lambda tva : State(np.array([tva[0]+tva[1], tva[0]-tva[1]])) # tva: Two-Valued Array
        forecast_state = lambda tva : tva + State(input=np.random.random(2) - 0.5)
        f = Filter(
            initial_state = [0,0],
            forecast_obs = forecast_obs,
            forecast_state = forecast_state,
        )
        f.load_measurements(50*[[2,0]])
        f.tools.set_covariance(0.01, 'measure')
        f.filter()
        print(f.get_last_state().round(1))
        assert f.get_last_state().round(1) == State([1,1])

    def test_noisy_tracking(self):
        """ Test: noisy sine

        Noisy data is filtered using sine's order 3 Taylor expansion
        """

        dt = 0.2

        def forecast_state(state):
            """ Using the raw approximation sin(x+h)~=sin(x)+sin(h) """
            return state + np.sin(dt)

        real_trajectory = np.sin(np.arange(0,10,dt))
        data = list(real_trajectory + np.random.normal(0,0.2,len(real_trajectory)))

        f = Filter(
                initial_state = 0,
                forecast_state = forecast_state
                )
        f.load_measurements(data)
        f.tools.set_covariance(0.2, 'measure')
        f.filter()

        rerr_ar = np.linalg.norm(f.extract_analyzed_states() \
            - real_trajectory)/np.linalg.norm(f.extract_analyzed_states()) #relative error, analyzed vs real
        rerr_dr = np.linalg.norm(np.array(data) \
            - real_trajectory)/np.linalg.norm(np.array(data)) #relative error, data vs real
        
        if True:
            import matplotlib.pyplot as plt 
            fig = plt.figure()
            plt.plot(np.array(data))
            plt.plot(f.extract_analyzed_states())
            plt.plot(real_trajectory)
            plt.legend(['data','analyzed','real'])
            plt.show()
            fig.savefig("tracking.pdf", bbox_inches='tight')
            print("relative error, analyzed vs real : " + str(rerr_ar))
            print("relative error, data vs     real : " + str(rerr_dr))

        assert rerr_ar < rerr_dr

    def test_constant_velocity_estimation(self):
        """ Test: constant velocity particle

        Noisy position is measured
        Trajectory is filtered and velocity is estimated
        """

        N = 100
        unc = 0.1
        exact_velocity = np.random.random()
        data = exact_velocity*np.arange(N) + np.random.normal(0,np.sqrt(exact_velocity),N)

        def forecast_state(state): # velocity [1] is constant, next position [0] is [0]+[1]
            return state @ np.array([[1,1],[0,1]])
        def forecast_obs(state):
            return State(state.get_values()[0])

        f = Filter(
                initial_state = [0.,1.], # initial pos: 0, initial speed: 1
                forecast_state = forecast_state,
                forecast_obs = forecast_obs
                )
        f.load_measurements(data)
        f.tools.set_covariance(unc, 'measure')
        f.filter()
        
        if True:
            import matplotlib.pyplot as plt 
            fig, axs = plt.subplots(2)
            #fig = plt.figure()
            axs[0].plot(f.extract_analyzed_states()[:,0])
            axs[0].plot(data)
            axs[0].plot(np.arange(N)*exact_velocity)
            axs[0].legend(['analyzed', 'measured', 'real'])
            #fig.savefig("trajectory.pdf", bbox_inches='tight')
            #figg = plt.figure()
            axs[1].plot(f.extract_analyzed_states()[:,1])
            axs[1].plot(exact_velocity*np.ones(N))
            axs[1].legend(['velocity estimate', 'real velocity'])
            #figg.savefig("velocity.pdf", bbox_inches='tight')
            plt.show()
        
        rerr = abs(exact_velocity - f.get_last_state().get_values()[1])/exact_velocity
        print(rerr)
        assert rerr < 0.1