# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: NIST regresion of cp

import numpy as np
import matplotlib.pyplot as plt

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
	plt.legend(["NIST", "fitted"])
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

	print(coeff_dH[7])
	print(dH_coeff)
	print(dH_diff)

	plt.figure(2)
	plt.plot(T, dH_calc, T, dH_fit)
	plt.legend(["NIST", "fitted"])
	plt.show()

	return dH_coeff

