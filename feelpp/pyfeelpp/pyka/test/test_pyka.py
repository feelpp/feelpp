from pyka import *
import numpy as np
import pytest

class TestState:

    def test_init(self):
        print('Testing State initialization...')
        s = State(input = [1,2,3,4,5])
        assert s.get_dim() == 5
        assert np.array_equal(s.get_values(), np.array([1,2,3,4,5]))
        assert str(s) == "State of dimension 5"

    def test_set_dim(self):
        s = State(dim = 3)
        assert np.array_equal(s.get_values(), np.array(3*[0]))

    def test_addition(self):
        print('Testing state addition...')
        assert State(input = [1,2]) + State(input = [2,3]) == State(input = [3,5])

    def test_multiplication(self):
        print('Testing scalar multiplication...')
        assert 3 * State(input = [1,2,3]) == State(input = [3,6,9])

    def test_weighted_sum(self):
        print('Testing state weighted sum')
        list = [State(input = [1,2]), State(input = [2,3]), State(input = [3,-4])]
        weights = [2, -1, 1]
        assert weighted_sum(list, weights) == State(input = [3, -3])

class TestFilter():

    def test_init(self):
        f = Filter(initial_state = [1,1,1])
        assert f.get_last_state().get_dim() == 3
        assert np.array_equal(f.get_last_state().get_values(), np.array([1,1,1]))
        assert str(f) == "Filter at time step 0\n    State dimension = 3\n    no observations"

    def test_load_data(self):
        f = Filter()
        data = [np.array([1,2,3]),np.array([2,3,4]),np.array([3,4,5]),np.array([4,5,6]),np.array([5,6,7])]
        f.load_measurements(data)
        assert f.real_observations[0] == State(input = [1,2,3])
        assert f.max_ts == 5

    def test_simple_dircet_filter(self):
        f = Filter(initial_state = 0)
        f.load_measurements(100*[1])
        f.filter()
        assert f.get_last_state().round(10) == State(input = [1.]).round(10)

    def test_multiple_direct_filter(self):
        dim = 1 + np.random.randint(100)
        start = np.random.random(dim)
        goal = np.random.random(dim)
        f = Filter(initial_state = start)
        f.load_measurements(100*[goal])
        f.filter()
        assert f.get_last_state().round(10) == State(input = goal).round(10)

    def test_sum_diff(self):
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
        print(f.get_last_state().round(3).get_values())
        assert f.get_last_state().round(1) == State([1,1])
