from python_things import ode_funcs
import numpy as np

def test_derivative():
    N_SPECIES = 100
    EXP_INTRA = 1
    EXP_INTER = 1
    K = 1
    a = np.loadtxt("a.txt")
    N_INITIAL = np.repeat(0.1,N_SPECIES)
    R_CONST = np.repeat(1,N_SPECIES)
    d = ode_funcs.derivative(0,N_INITIAL,R_CONST, a,K,EXP_INTRA,EXP_INTER)
    print(d)

test_derivative()
    