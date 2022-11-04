from ._timing import *

class Timer(object):
    """helper class to time a block of code

    Args:
        object (_type_): object to time
    """
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print('[%s]' % self.name)
        elapsed=time.time() - self.tstart
        print(f'Elapsed: {elapsed}]')
        return  elapsed