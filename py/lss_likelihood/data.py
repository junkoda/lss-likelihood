import lss_likelihood._lss_likelihood as c
import numpy as np

class Data:
    def __init__(self, P0, P2, P4):
        """
        Args:
          P0, P2, P4 (array)
        """

        self._data = c._data_alloc(P0, P2, P4)

