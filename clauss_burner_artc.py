# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: Calculate burner with those equations:
#			   H2S + 3/2 O2  ->    SO2 +   H2O
#			 2 H2S +    SO2 <=> 3/2 S2 + 2 H2O
#			   Positions in vector: H2S O2 SO2 H2O S2 CH4 CS2 CO C2H6 C3H8 CO2 COS NH3 N2 H2 Ar T
#			   Positions in vector: H2S, SO2, S2, H2O, O2, N2, T

import numpy as np 

def clauss_burner_artc(yvec, V, pars):
	R = 8.314
	p = pars

	n    = np.empty(len(yvec-1))
	y    = np.empty(len(yvec))
	c    = np.empty(len(yvec))
	cp   = np.empty(len(yvec))
	r    = np.empty(2)
	drH  = np.empty(len(r))
	k    = np.empty(len(r))
	K    = np.empty(len(r))

	n[0]  = yvec[0]
	n[1]  = yvec[1]
	n[2]  = yvec[2]
	n[3]  = yvec[3]
	n[4]  = yvec[4]
	n[5]  = yvec[5]
	n[6]  = yvec[6]
	n[7]  = yvec[7]
	n[8]  = yvec[8]
	n[9]  = yvec[9]
	n[10] = yvec[10]
	n[11] = yvec[11]
	n[12] = yvec[12]
	n[13] = yvec[13]
	n[14] = yvec[14]
	n[15] = yvec[15]
	T     = yvec[16]

	nTOT = n[0] + n[1] + n[2] + n[3] + n[4] + n[5] + n[6] + n[7] + n[8] + n[9] + n[10] + n[11] + n[12] + n[13] + n[14] + n[15]

	y[0]  = n[0]  / nTOT
	y[1]  = n[1]  / nTOT
	y[2]  = n[2]  / nTOT
	y[3]  = n[3]  / nTOT
	y[4]  = n[4]  / nTOT
	y[5]  = n[5]  / nTOT
	y[6]  = n[6]  / nTOT
	y[7]  = n[7]  / nTOT
	y[8]  = n[8]  / nTOT
	y[9]  = n[9]  / nTOT
	y[10] = n[10] / nTOT
	y[11] = n[11] / nTOT
	y[12] = n[12] / nTOT
	y[13] = n[13] / nTOT
	y[14] = n[14] / nTOT
	y[15] = n[15] / nTOT

	c[0]  = y[0]  * p/(R*T)
	c[1]  = y[1]  * p/(R*T)
	c[2]  = y[2]  * p/(R*T)
	c[3]  = y[3]  * p/(R*T)
	c[4]  = y[4]  * p/(R*T)
	c[5]  = y[5]  * p/(R*T)
	c[6]  = y[6]  * p/(R*T)
	c[7]  = y[7]  * p/(R*T)
	c[8]  = y[8]  * p/(R*T)
	c[9]  = y[9]  * p/(R*T)
	c[10] = y[10] * p/(R*T)
	c[11] = y[11] * p/(R*T)
	c[12] = y[12] * p/(R*T)
	c[13] = y[13] * p/(R*T)
	c[14] = y[14] * p/(R*T)
	c[15] = y[15] * p/(R*T)

	cp[0] = 31.94  + 0.001436*T + 0.00002432*pow(T, 2) - 1.176e-8*pow(T, 3)
	cp[4] = 47.7   + 0.007171*T - 0.00008556/pow(T, 2)
	cp[1] = 36.162 + 0.000845*T - 0.0000431 /pow(T, 2)
	cp[3] = 30.12  + 0.0113  *T
	cp[2] = 35.0
	cp[5] = 0.00418 *T                                             + 27.82

	drH[0] = -(2.31198e-35*(2.19016e40*T + 1.60936e36*pow(T, 2) + 1.02726e30*pow(T, 3) + 3.50637e29*pow(T, 4) - 1.27164e26*pow(T, 5) + 3.68688e30))/T
	drH[1] =  (4.62396e-35*(9.48365e38*T + 2.74613e35*pow(T, 2) + 2.04186e32*pow(T, 3) - 3.50637e29*pow(T, 4) + 1.27164e26*pow(T, 5) - 9.32101e29))/T

	k[0] = 1*np.exp(1000*(1/500 - 1/T))
	k[1] = 1*np.exp(1000*(1/500 - 1/T))

	K[0] = np.exp(1.17873e-10*pow(T, 3) - 4.47534*np.log(T) - 4.87531e-7*pow(T, 2) - 0.00000285663*T + (60904.6*T + 0.00000512629)/pow(T, 2) + 21.1096)
	K[1] = np.exp(0.00113561*T + 1.5273*np.log(T) - 9.75062e-7*pow(T, 2) + 2.35747e-10*pow(T, 3) - (1.0*(5274.48*T - 0.00000259201))/pow(T, 2) - 3.1067)

	r[0] = k[0]*(    c[0]   * pow(c[4], 3/2)  -     c[1]      *     c[3]     / K[0])
	r[1] = k[1]*(pow(c[0],2)*     c[1]        - pow(c[2],3/2) * pow(c[3], 2) / K[1])

	dnH2SdV = - 1.0 * r[0] - 2.0 * r[1]
	dnSO2dV = + 1.0 * r[0] - 1.0 * r[1]
	dnS2dV  =              + 3/2 * r[1]
	dnH2OdV = + 1.0 * r[0] + 2.0 * r[1]
	dnO2dV  = - 3/2 * r[0]
	dnN2dV  =   0

	dTdV = (- r[0]*drH[0] - r[1]*drH[1]) / (n[0]*cp[0] + n[1]*cp[1] + n[2]*cp[2] + n[3]*cp[3] + n[4]*cp[4] + n[5]*cp[5])

	dydt = np.empty(len(yvec))
	
	dydt[0] = dnH2SdV
	dydt[1] = dnSO2dV
	dydt[2] = dnS2dV
	dydt[3] = dnH2OdV
	dydt[4] = dnO2dV
	dydt[5] = dnN2dV
	dydt[6] = dTdV

	return dydt 
