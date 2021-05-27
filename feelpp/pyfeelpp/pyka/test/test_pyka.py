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
        f = Filter(initial_state = 0)
        f.load_measurements(100*[1])
        f.filter()
        assert f.get_last_state().round(10) == State(input = [1.]).round(10)

    def test_multiple_direct_filter(self):
        """ Test: direct nD filter with identity forecasts """
        dim = 1 + np.random.randint(100)
        start = np.random.random(dim)
        goal = np.random.random(dim)
        f = Filter(initial_state = start)
        f.load_measurements(100*[goal])
        f.filter()
        assert f.get_last_state().round(10) == State(input = goal).round(10)

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
        assert f.get_last_state().round(1) == State([1,1])

    def test_noisy_tracking(self):

        dt = 0.05

        def forecast_state(state):
            """ Using the raw approximation sin(x+h)~=sin(x)+sin(h) """
            return state + np.sin(dt)

        real_trajectory = np.sin(np.arange(0,10,0.1))
        data = list(real_trajectory + np.random.normal(0,0.2,len(real_trajectory)))

        f = Filter(
                initial_state = 0,
                forecast_state = forecast_state
                )
        f.load_measurements(data)
        f.tools.set_covariance(0.5, 'measure')
        f.filter()

        rerr_ad = np.linalg.norm(f.extract_analyzed_states() \
            - np.array(data))/np.linalg.norm(f.extract_analyzed_states()) #relative error, analyzed vs data
        rerr_rd = np.linalg.norm(f.extract_analyzed_states() \
            - real_trajectory)/np.linalg.norm(f.extract_analyzed_states()) #relative error, real vs data
        
        if False:
            import matplotlib.pyplot as plt 
            plt.plot(np.array(data))
            plt.plot(f.extract_analyzed_states())
            plt.plot(real_trajectory)
            plt.show()
            print("relative error, analyzed vs data : " + str(np.linalg.norm(f.extract_analyzed_states() - np.array(data))/np.linalg.norm(f.extract_analyzed_states())))
            print("relative error, real vs     data : " + str(np.linalg.norm(f.extract_analyzed_states() - real_trajectory)/np.linalg.norm(f.extract_analyzed_states())))

        assert rerr_ad < rerr_rd
