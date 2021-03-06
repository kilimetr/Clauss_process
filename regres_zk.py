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

# from fit_regres import poly_regres_cp, poly_regres_dH

Tmin  = 298
Tmax  = 1400
Tstep = round((1400-298)/25)

Trange = [Tmin, Tmax, Tstep]

order = 3

def poly_regres_cp(coeff_cp, Trange, order):
	T = np.linspace(Trange[0], Trange[1], Trange[2])

	TT = T/1000 # NIST condition

	cp_calc = coeff_cp[0] + coeff_cp[1]*pow(TT, 1) + coeff_cp[2]*pow(TT, 2) + coeff_cp[3]*pow(TT, 3) + coeff_cp[4]/pow(TT, 2)

	cp_coeff = np.polyfit(T, cp_calc, order)

	cp_fit = cp_coeff[3] + cp_coeff[2]*pow(T, 1) + cp_coeff[1]*pow(T, 2) + cp_coeff[0]*pow(T, 3)


	cp_diff = np.absolute(cp_calc - cp_fit)
	cp_diff = sum(cp_diff)

	print(cp_coeff)
	print(cp_diff)

	plt.figure(1)
	plt.plot(T, cp_calc, T, cp_fit)
	plt.show()

	return cp_coeff


def poly_regres_dH(coeff_dH, Trange, order):
	T = np.linspace(Trange[0], Trange[1], Trange[2])

	TT = T/1000

	dH_calc = coeff_dH[0]*TT + coeff_dH[1]*pow(TT, 2)/2 + coeff_dH[2]*pow(TT, 3)/3 + coeff_dH[3]*pow(TT, 4)/4 - coeff_dH[4]/TT + coeff_dH[5] - coeff_dH[6]
	
	dH_coeff = np.polyfit(T, dH_calc, order)

	dH_fit = dH_coeff[3] + dH_coeff[2]*T + dH_coeff[1]*pow(T, 2) + dH_coeff[0]*pow(T, 3)

	dH_diff = np.absolute(dH_fit - dH_calc)
	dH_diff = sum(dH_diff)

	print(dH_coeff)
	print(dH_diff)

	plt.figure(2)
	plt.plot(T, dH_calc, T, dH_fit)
	plt.show()

	return dH_coeff


A  = 26.88412
B  = 18.67809
C  = 3.434203
D  = -3.378702	
E  = 0.135882
F  = -28.91211
H  = -20.50202
dH = -20.50
coef_cp = [A, B, C, D, E]
coef_dH = [A, B, C, D, E, F, H]

H2S_coef_cp = poly_regres_cp(coef_cp, Trange, order)

H2S_coef_dH = poly_regres_dH(coef_dH, Trange, order)

cancel

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




def poly_regres_cp(coeff_cp, Trange, order):
	T = np.linspace(Trange[0], Trange[1], Trange[2])

	TT = T/1000 # NIST condition

	cp_calc = coeff_cp[0] + coeff_cp[1]*pow(TT, 1) + coeff_cp[2]*pow(TT, 2) + coeff_cp[3]*pow(TT, 4) + coeff_cp[4]*pow(TT, 5)

	cp_coeff = np.polyfit(T, cp_calc, order)

	cp_fit = cp_coeff[3] + cp_coeff[2]*pow(T, 1) + cp_coeff[1]*pow(T, 2) + cp_coeff[0]*pow(T, 3)

	cp_diff = np.absolute(cp_calc - cp_fit)

	print(cp_coeff)
	print(cp_diff)

	plt.figure(1)
	plt.subplot(2,1,1)
	plt.plot(T, cp_calc)
	plt.subplot(2,1,2)
	plt.plot(T, cp_fit)
	plt.show()

	return cp_coeff


def poly_regres_dH(coeff_dH, Trange, order):
	T = np.linspace(Trange[0], Trange[1], Trange[2])

	TT = T/1000

	dH_calc = coeff_dH[0]*TT + coeff_dH[1]*pow(TT, 2)/2 + coeff_dH[2]*pow(TT, 3)/3 + coeff_dH[3]*pow(TT, 4) - coeff_dH[4]/TT + coeff_dH[5] - coeff_dH[6] + coeff_dH[7]
	dH_calc = dH_calc * 1000 # [J/mol]
	
	dH_coeff = np.polyfit(T, dH_calc, order)

	dH_fit = dH_coeff[3] + dH_coeff[2]*T + dH_coeff[1]*pow(T, 2) + dH_coeff[0]*pow(T, 3)

	dH_diff = np.absolute(dH_fit - dH_calc)

	print(dH_coeff)
	print(dH_diff)

	plt.figure(2)
	plt.subplot(2, 1, 1)
	plt.plot(T, dH_calc)
	plt.subplot(2, 1, 2)
	plt.plot(T, dH_fit)
	plt.show()

	return dH_coeff


