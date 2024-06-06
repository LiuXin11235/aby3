import numpy as np
import math

def get_k(V, E, B=1024):
    """compute the proper k value, supposing the edges are distributed evenly among the edge blocks.
    """
    k = int(np.ceil(B*V/E))
    k = 2**math.floor(math.log2(k))

    return k
