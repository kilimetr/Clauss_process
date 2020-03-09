# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: NIST regresion of cp

import numpy as np
import matplotlib.pyplot as plt


def poly_regres_cp(coeff, Trange, order):
	T = np.linspace(Trange[0], Trange[1], Trange[2])

	TT = T/1000 # NIST condition

	value_calc = coeff[0] + coeff[1]*pow(TT, 1) + coeff[2]*pow(TT, 2) + coeff[3]*pow(TT, 4) + coeff[4]*pow(TT, 5)

	value_coeff = np.polyfit(T, value_calc, order)

	value_fit = value_coeff[3] + value_coeff[2]*pow(T, 1) + value_coeff[1]*pow(T, 2) + value_coeff[0]*pow(T, 3)

	value_diff = np.absolute(value_calc - value_fit)

	print(value_coeff)
	print(value_diff)

	plt.figure(1)
	plt.subplot(2,1,1)
	plt.plot(T, value_calc)
	plt.subplot(2,1,2)
	plt.plot(T, value_fit)
	plt.show()

	return value_coeff


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

