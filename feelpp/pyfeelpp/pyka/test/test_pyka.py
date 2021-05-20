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
        f.load_real_observations(data)
        assert f.real_observations[0] == State(input = [1,2,3])
        assert f.max_ts == 5

    def test_simple_filter(self):
        f = Filter()
        f.filter(initial_guess = 0, data = 100*[1])
        assert f.get_last_state() == State(input = [1.])
