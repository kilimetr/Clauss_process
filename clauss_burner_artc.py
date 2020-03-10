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

	cp[0]  = 31.94 + 0.001436*T + 0.00002432*pow(T, 2) - 1.176e-08*pow(T, 3)
	cp[1]  = 47.7  + 0.007171*T - 0.00008556/pow(T, 2)
	cp[2]  = 36.16 + 0.000845*T - 0.0000431 /pow(T, 2)
	cp[3]  = 3.01308973e+01 + 8.36012484e-03*T + 1.02417826e-06*pow(T, 2) - 5.99505721e-11*pow(T, 3)
	cp[4]  = 35.0

	cp[5]  = 2.10774412e+01 + 1.85050281e-02*T + 1.21547535e-05*pow(T, 2) + 8.45010505e-09*pow(T, 3)
	cp[6]  = 4.33539314e+01 + 1.42933993e-02*T + 3.52004019e-06*pow(T, 2) + 1.43302897e-09*pow(T, 3)
	cp[7]  = 2.59045855e+01 + 6.17807249e-03*T - 3.17133389e-08*pow(T, 2) - 7.28804158e-10*pow(T, 3)

	cp[8]  = # C2H6
	cp[9]  = # C3H8
	cp[10] = 3.42873947e+01 + 1.27694188e-02*T + 4.10801203e-06*pow(T, 2) + 2.11639058e-09*pow(T, 3)
	cp[11] = 4.10871182e+01 + 1.25318691e-02*T + 2.41945066e-06*pow(T, 2) + 6.77642806e-10*pow(T, 3)
	cp[12] = 2.89461375e+01 + 1.32137091e-02*T + 5.66727810e-06*pow(T, 2) + 3.41943966e-09*pow(T, 3)
	cp[13] = 2.22747069e+01 + 7.43818001e-03*T + 1.89387351e-06*pow(T, 2) + 7.99616401e-10*pow(T, 3)
	cp[14] = 2.93021725e+01 + 5.25410813e-03*T - 1.60316925e-06*pow(T, 2) - 1.98798627e-09*pow(T, 3)
	cp[15] = 1.95675260e+01 + 4.26884125e-03*T - 3.83532582e-07*pow(T, 2) - 8.17473952e-10*pow(T, 3)

	/Users/kilimetr/Desktop/python/Clauss_process/regres.py:37: RankWarning: Polyfit may be poorly conditioned
  H2S_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[5.11477506e-10 1.78306247e-06 9.12923232e-03 2.98364463e+01]
[3.55271368e-14 1.42108547e-14]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:38: RankWarning: Polyfit may be poorly conditioned
  H2S_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 6.70367713e-06  8.20354966e-03  3.67545353e+00 -2.25228053e+04]
[1.45519152e-11 3.63797881e-12]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:53: RankWarning: Polyfit may be poorly conditioned
  O2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[-1.49534151e-09 -8.94616113e-07  6.33239164e-03  3.05289723e+01]
[2.84217094e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:54: RankWarning: Polyfit may be poorly conditioned
  O2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 5.78093859e-06  7.68747370e-03  7.85808005e+00 -2.70324897e+03]
[2.21689334e-12 1.45519152e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:69: RankWarning: Polyfit may be poorly conditioned
  SO2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[3.53971581e-09 5.99930802e-06 1.46931655e-02 3.35826314e+01]
[4.26325641e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:70: RankWarning: Polyfit may be poorly conditioned
  SO2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 3.09962707e-05  3.11201895e-02 -3.50909552e+01 -2.89877255e+05]
[2.91038305e-10 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:85: RankWarning: Polyfit may be poorly conditioned
  H2O_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[-5.99505721e-11  1.02417826e-06  8.36012484e-03  3.01308973e+01]
[2.84217094e-14 7.10542736e-15]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:86: RankWarning: Polyfit may be poorly conditioned
  H2O_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 1.80135000e-05  1.57441057e-02 -3.82980803e+01 -2.32311516e+05]
[8.73114914e-11 1.74622983e-10]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:101: RankWarning: Polyfit may be poorly conditioned
  CH4_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[8.45010505e-09 1.21547535e-05 1.85050281e-02 2.10774412e+01]
[4.61852778e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:102: RankWarning: Polyfit may be poorly conditioned
  CH4_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 1.68489876e-05  1.98339154e-02  3.23609899e+00 -7.80110469e+04]
[7.27595761e-11 3.09228199e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:117: RankWarning: Polyfit may be poorly conditioned
  CS2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[1.43302897e-09 3.52004019e-06 1.42933993e-02 4.33539314e+01]
[4.97379915e-14 1.42108547e-14]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:118: RankWarning: Polyfit may be poorly conditioned
  CS2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[8.44749258e-06 1.51579458e-02 4.14937008e+01 1.03069097e+05]
[1.30967237e-10 2.91038305e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:133: RankWarning: Polyfit may be poorly conditioned
  CO_coef_cp = poly_regres_cp(coeff, Trange, order)
[-7.28804158e-10 -3.17133389e-08  6.17807249e-03  2.59045855e+01]
[2.13162821e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:134: RankWarning: Polyfit may be poorly conditioned
  CO_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 9.89762065e-06  9.40442621e-03 -1.52792091e+01 -1.07090818e+05]
[1.16415322e-10 4.36557457e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:149: RankWarning: Polyfit may be poorly conditioned
  CO2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[2.11639058e-09 4.10801203e-06 1.27694188e-02 3.42873947e+01]
[4.97379915e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:150: RankWarning: Polyfit may be poorly conditioned
  CO2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 3.22828525e-05  2.95636828e-02 -5.83281621e+01 -3.79579379e+05]
[3.49245965e-10 5.82076609e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:165: RankWarning: Polyfit may be poorly conditioned
  COS_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[6.77642806e-10 2.41945066e-06 1.25318691e-02 4.10871182e+01]
[4.97379915e-14 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:166: RankWarning: Polyfit may be poorly conditioned
  COS_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 1.88385650e-05  2.03300454e-02 -1.04976462e+01 -1.37560122e+05]
[5.82076609e-11 3.63797881e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:181: RankWarning: Polyfit may be poorly conditioned
  NH3_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[3.41943966e-09 5.66727810e-06 1.32137091e-02 2.89461375e+01]
[7.10542736e-15 2.84217094e-14]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:182: RankWarning: Polyfit may be poorly conditioned
  NH3_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 1.17575494e-05  1.40540446e-02  3.89128386e+00 -4.86140729e+04]
[4.36557457e-11 1.81898940e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:196: RankWarning: Polyfit may be poorly conditioned
  N2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[7.99616401e-10 1.89387351e-06 7.43818001e-03 2.22747069e+01]
[2.13162821e-14 7.10542736e-15]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:197: RankWarning: Polyfit may be poorly conditioned
  N2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 8.22929145e-06  9.36925322e-03 -8.50606335e-01 -4.67712832e+04]
[4.36557457e-11 1.18234311e-11]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:212: RankWarning: Polyfit may be poorly conditioned
  H2_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[-1.98798627e-09 -1.60316925e-06  5.25410813e-03  2.93021725e+01]
[4.26325641e-14 1.42108547e-14]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:213: RankWarning: Polyfit may be poorly conditioned
  H2_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 3.76336698e-06  4.99340479e-03  5.03062895e+00 -2.06273070e+03]
[3.62376795e-12 0.00000000e+00]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:228: RankWarning: Polyfit may be poorly conditioned
  Ar_coef_cp = poly_regres_cp(coef_cp, Trange, order)
[-8.17473952e-10 -3.83532582e-07  4.26884125e-03  1.95675260e+01]
[1.77635684e-14 3.55271368e-15]
/Users/kilimetr/Desktop/python/Clauss_process/regres.py:229: RankWarning: Polyfit may be poorly conditioned
  Ar_coef_dH = poly_regres_dH(coef_dH, Trange, order)
[ 3.70590196e-06  4.91785135e-03  4.95911842e+00 -2.01231587e+03]
[1.22923893e-12 3.63797881e-12]
[Finished in 10.6s]


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
