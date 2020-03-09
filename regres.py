# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: NIST component regresion 


import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/kilimetr/Desktop/python/Clauss_process")

from fit_regres import poly_regres_cp

Tmin  = 298
Tmax  = 1400
Tstep = 2

Trange = [Tmin, Tmax, Tstep]

order = 3


A = 30.09200
B = 6.832514
C = 6.793435
D = -2.534480
E = 0.082139
coeff = [A, B, C, D, E]

H2O_coef = poly_regres_cp(coeff, Trange, order)


A = 19.50583
B = 19.88705
C = -8.598535
D = 1.369784
E = 0.527601
coeff = [A, B, C, D, E]

N2_coef = poly_regres_cp(coeff, Trange, order)


A = -0.703029
B = 108.4773
C = -42.52157
D = 5.862788
E = 0.678565
coeff = [A, B, C, D, E]

CH4_coef = poly_regres_cp(coeff, Trange, order)


A = 35.85391
B = 52.49121
C = -40.83743
D = 12.00155
E = -0.224831
coeff = [A, B, C, D, E]

CS2_coef = poly_regres_cp(coeff, Trange, order)


A = 25.56759
B = 6.096130
C = 4.054656
D = -2.671301
E = 0.131021
coeff = [A, B, C, D, E]

CO_coef = poly_regres_cp(coeff, Trange, order)


A = 24.99735
B = 55.18696
C = -33.69137
D = 7.948387
E = -0.136638
coeff = [A, B, C, D, E]

CO2_coef = poly_regres_cp(coeff, Trange, order)


A = 34.53892
B = 43.05378
C = -26.61773
D = 6.338844
E = -0.327515
coeff = [A, B, C, D, E]

COS_coef = poly_regres_cp(coeff, Trange, order)


A = 19.99563
B = 49.77119
C = -15.37599
D = 1.921168
E = 0.189174
coeff = [A, B, C, D, E]

NH3_coef = poly_regres_cp(coeff, Trange, order)


A = 33.066178
B = -11.363417
C = 11.432816
D = -2.772874
E = -0.158558
coeff = [A, B, C, D, E]

H2_coef = poly_regres_cp(coeff, Trange, order)


A = 20.78600
B = pow(2.825911, -7)
C = -pow(1.464191, -7)
D = pow(1.092131, -8)
E = -pow(3.661371, -8)
coeff = [A, B, C, D, E]

Ar_coef = poly_regres_cp(coeff, Trange, order)



