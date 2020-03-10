# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: NIST regresion of cp

import numpy as np
import matplotlib.pyplot as plt


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

# import numpy as np
# x = np.array([0,1,2,3,4,5])
# y = np.array([0,0.8,0.9,0.1,-0.8,-1])
# print(x)
# print(y)

# p1 = np.polyfit(x,y,1)
# p2 = np.polyfit(x,y,2)
# p3 = np.polyfit(x,y,3)
# print(p1)
# print(p2)
# print(p3)

# import matplotlib.pyplot as plt
# plt.plot(x,y,'o')
# xp = np.linspace(-2,6,100)
# plt.plot(xp,np.polyval(p1,xp),'r-')
# plt.plot(xp,np.polyval(p2,xp),'b--')
# plt.plot(xp,np.polyval(p3,xp),'m:')
# yfit = p1[0] * x + p1[1]
# yresid= y - yfit
# SSresid = np.sum(yresid**2)
# SStotal = len(y) * np.var(y)
# rsq = 1 - SSresid/SStotal
# print(yfit)
# print(y)
# print(rsq)

# from scipy.stats import linregress
# slope,intercept,r_value,p_value,std_err = linregress(x,y)
# print(r_value**2)
# print(p_value)
# plt.show()

