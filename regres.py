# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: NIST cp and dH component regresion
# thermal capacity    cp    [J/mol/K]
# Enthalpy            dH    [J/mol]


import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir("/Users/kilimetr/Desktop/python/Clauss_process")

from fit_regres import poly_regres_cp, poly_regres_dH

Tmin  = 298
Tmax  = 1400
Tstep = 2

Trange = [Tmin, Tmax, Tstep]

order = 3


A  = 26.88412
B  = 18.67809
C  = 3.434203
D  = -3.378702	
E  = 0.135882
F  = -28.91211
H  = -20.50202
dH = -20.50
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

H2S_coef_cp = poly_regres_cp(coef_cp, Trange, order)
H2S_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 30.03235
B  = 8.772972
C  = -3.988133
D  = 0.788313
E  = -0.741599
F  = -11.32468
H  = 0.0
dH = 0
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

O2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
O2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A = 21.43049
B = 74.35094
C = -57.75217
D = 16.35534
E = 0.086731
F = -305.7688
H = -296.8422
dH = -296.84
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

SO2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
SO2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 30.09200
B  = 6.832514
C  = 6.793435
D  = -2.534480
E  = 0.082139
F  = -250.8810
H  = -241.8264
dH = -241.83
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

H2O_coef_cp = poly_regres_cp(coef_cp, Trange, order)
H2O_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = -0.703029
B  = 108.4773
C  = -42.52157
D  = 5.862788
E  = 0.678565
F  = -76.84376
H  = -74.87310
dH = -74.87
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

CH4_coef_cp = poly_regres_cp(coef_cp, Trange, order)
CH4_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 35.85391
B  = 52.49121
C  = -40.83743
D  = 12.00155
E  = -0.224831
F  = 103.5030
H  = 116.9432
dH = 116.94
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

CS2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
CS2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 25.56759
B  = 6.096130
C  = 4.054656
D  = -2.671301
E  = 0.131021
F  = -118.0089
H  = -110.5271
dH = -110.53
coeff = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

CO_coef_cp = poly_regres_cp(coeff, Trange, order)
CO_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 24.99735
B  = 55.18696
C  = -33.69137
D  = 7.948387
E  = -0.136638
F  = -403.6075
H  = -393.5224
dH = -393.52
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

CO2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
CO2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 34.53892
B  = 43.05378
C  = -26.61773
D  = 6.338844
E  = -0.327515
F  = -151.5001
H  = -138.4071
dH = -138.41
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

COS_coef_cp = poly_regres_cp(coef_cp, Trange, order)
COS_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 19.99563
B  = 49.77119
C  = -15.37599
D  = 1.921168
E  = 0.189174
F  = -53.30667
H  = -45.89806
dH = -45.90
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

NH3_coef_cp = poly_regres_cp(coef_cp, Trange, order)
NH3_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 19.50583
B  = 19.88705
C  = -8.598535
D  = 1.369784
E  = 0.527601
F  = -4.935202
H  = 0
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

N2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
N2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 33.066178
B  = -11.363417
C  = 11.432816
D  = -2.772874
E  = -0.158558
F  = -9.980797
H  = 0
dH = 0
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

H2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
H2_coef_dH = poly_regres_dH(coef_dH, Trange, order)



A  = 20.78600
B  = pow(2.825911, -7)
C  = -pow(1.464191, -7)
D  = pow(1.092131, -8)
E  = -pow(3.661371, -8)
F  = -6.197350
H  = 0
dH = 0
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, dH]

Ar_coef_cp = poly_regres_cp(coef_cp, Trange, order)
Ar_coef_dH = poly_regres_dH(coef_dH, Trange, order)


