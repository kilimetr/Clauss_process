
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/kilimetr/Desktop/python/Clauss_process")

from thermo import poly_regres

A = 30.09200
B = 6.832514
C = 6.793435
D = -2.534480
E = 0.082139

Tmin  = 298
Tmax  = 1400
Tstep = 2

coeff = [A, B, C, D, E]
Trange = [Tmin, Tmax, Tstep]

order = 3

coef = poly_regres(coeff, Trange, order)