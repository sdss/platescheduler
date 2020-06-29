

class Cadence(object):
    """Base class for cadences
    """

    def __init__(self, delta=0):
        self.delta = delta


    def __call__(self, since_previous):
        return self.delta > since_previous
