import numpy as np
import numpy.polynomial.hermite_e as He

def py_hermite_e(x, n):
    
    x = np.asarray(x)
    scalar_input = False
    if x.ndim == 0:
        x = x[np.newaxis]  # Makes x 1D
        scalar_input = True

    c = [0 if d != n else 1 for d in range(0, n+1)]
    H = He.HermiteE(c)
    ret = [H(i) for i in x]
        
    if scalar_input:
        return np.squeeze(ret)
    return ret
