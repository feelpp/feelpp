from base import *

filt = Filter()
filt.set_state([4])
filt.load_real_observations([State([2])])

filt.forecast()