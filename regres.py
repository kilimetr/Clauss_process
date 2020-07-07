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

Tmin  = 400
Tmax  = 1400
Tstep = round((Tmax - Tmin)/25)

Trange = [Tmin, Tmax, Tstep]

order = 3

plottt = False


A  = 26.88412
B  = 18.67809
C  = 3.434203
D  = -3.378702	
E  = 0.135882
F  = -28.91211
H  = -20.50202
dH = -20.5
xxx = "H2S"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

H2S_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
H2S_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 30.03235
B  = 8.772972
C  = -3.988133
D  = 0.788313
E  = -0.741599
F  = -11.32468
H  = 0.0
dH = 0
xxx = "O2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

O2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
O2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A = 21.43049
B = 74.35094
C = -57.75217
D = 16.35534
E = 0.086731
F = -305.7688
H = -296.8422
dH = -296.84
xxx = "SO2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

SO2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
SO2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 30.09200
B  = 6.832514
C  = 6.793435
D  = -2.534480
E  = 0.082139
F  = -250.8810
H  = -241.8264
dH = -241.83
xxx = "H2O"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

H2O_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
H2O_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 33.51313
B  = 5.065360
C  = -1.059670
D  = 0.089905	
E  = -0.211911
F  = 117.6855
H  = 128.6003
dH = 128.60
xxx = "S2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

S2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
S2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = -0.703029
B  = 108.4773
C  = -42.52157
D  = 5.862788
E  = 0.678565
F  = -76.84376
H  = -74.87310
dH = -74.87
xxx = "CH4"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

CH4_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
CH4_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 35.85391
B  = 52.49121
C  = -40.83743
D  = 12.00155
E  = -0.224831
F  = 103.5030
H  = 116.9432
dH = 116.94
xxx = "CS2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

CS2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
CS2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 25.56759
B  = 6.096130
C  = 4.054656
D  = -2.671301
E  = 0.131021
F  = -118.0089
H  = -110.5271
dH = -110.53
xxx = "CO"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

CO_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
CO_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 24.99735
B  = 55.18696
C  = -33.69137
D  = 7.948387
E  = -0.136638
F  = -403.6075
H  = -393.5224
dH = -393.52
xxx = "CO2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

CO2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
CO2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 34.53892
B  = 43.05378
C  = -26.61773
D  = 6.338844
E  = -0.327515
F  = -151.5001
H  = -138.4071
dH = -138.41
xxx = "COS"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

COS_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
COS_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 19.99563
B  = 49.77119
C  = -15.37599
D  = 1.921168
E  = 0.189174
F  = -53.30667
H  = -45.89806
dH = -45.90
xxx = "NH3"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

NH3_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
NH3_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 19.50583
B  = 19.88705
C  = -8.598535
D  = 1.369784
E  = 0.527601
F  = -4.935202
H  = 0
xxx = "N2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

N2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
N2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 33.066178
B  = -11.363417
C  = 11.432816
D  = -2.772874
E  = -0.158558
F  = -9.980797
H  = 0
dH = 0
xxx = "H2"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

H2_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
H2_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)



A  = 20.78600
B  = pow(2.825911, -7)
C  = -pow(1.464191, -7)
D  = pow(1.092131, -8)
E  = -pow(3.661371, -8)
F  = -6.197350
H  = 0
dH = 0
xxx = "Ar"
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H, xxx]

Ar_coef_cp = poly_regres_cp(coef_cp, Trange, order, plottt)
Ar_coef_dH = poly_regres_dH(coef_dH, Trange, order, plottt)

results = {
	"H2S": [H2S_coef_cp, H2S_coef_dH],
	"O2":  [O2_coef_cp,  O2_coef_dH],
	"SO2": [SO2_coef_cp, SO2_coef_dH],
	"H2O": [H2O_coef_cp, H2O_coef_dH],
	"S2":  [S2_coef_cp,  S2_coef_dH],
	"CH4": [CH4_coef_cp, CH4_coef_dH],
	"CS2": [CS2_coef_cp, CS2_coef_dH],
	"CO":  [CO_coef_cp,  CO_coef_dH],
	"CO2": [CO2_coef_cp, CO2_coef_dH],
	"COS": [COS_coef_cp, COS_coef_dH],
	"NH3": [NH3_coef_cp, NH3_coef_dH],
	"N2":  [N2_coef_cp,  N2_coef_dH],
	"H2":  [H2_coef_cp,  H2_coef_dH],
	"Ar":  [Ar_coef_cp,  Ar_coef_dH]
}

print(results)