# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: Comparing thermo between two models

import numpy as np
import matplotlib.pyplot as plt

T     = np.linspace(298.15, 1500, 800)
cp_1  = np.empty([len(T), 6])
drH_1 = np.empty([len(T), 2])
cp_2  = np.empty([len(T), 6])
dH_2  = np.empty([len(T), 6])
drH_2 = np.empty([len(T), 6])

#	FIRST MODEL'S DATA
#			   H2S + 3/2 O2  ->    SO2 +   H2O
#			 2 H2S +    SO2 <=> 3/2 S2 + 2 H2O
#			   Positions in vector: H2S, SO2, S2, H2O, O2, N2


cp_1[:,0] = 0.001436*T + 0.00002432*pow(T, 2) - 1.176e-8*pow(T, 3) + 31.94		# H2S
cp_1[:,1] = 0.000845*T - 0.0000431 /pow(T, 2)                      + 36.162		# SO2
cp_1[:,2] =                                                        + 35.0		# S2
cp_1[:,3] = 0.0113  *T                                             + 30.12		# H2O
cp_1[:,4] = 0.007171*T - 0.00008556/pow(T, 2)                      + 47.7		# O2
cp_1[:,5] = 0.00418 *T                                             + 27.82		# N2

drH_1[:,0] = -(2.31198e-35*(2.19016e40*T + 1.60936e36*pow(T, 2) + 1.02726e30*pow(T, 3) + 3.50637e29*pow(T, 4) - 1.27164e26*pow(T, 5) + 3.68688e30))/T
drH_1[:,1] =  (4.62396e-35*(9.48365e38*T + 2.74613e35*pow(T, 2) + 2.04186e32*pow(T, 3) - 3.50637e29*pow(T, 4) + 1.27164e26*pow(T, 5) - 9.32101e29))/T



#	SECOND MODEL'S DATA
#			   H2S + 3/2 O2  ->    SO2 +   H2O
#			 2 H2S +    SO2 <=> 3/2 S2 + 2 H2O
#			   Positions in vector: H2S O2 SO2 H2O S2 N2

cp_2[:,0]  = 2.98364463e+01 + 9.12923232e-03*T + 1.78306247e-06*pow(T, 2) + 5.11477506e-10*pow(T, 3)	# H2S
cp_2[:,1]  = 3.05289723e+01 + 6.33239164e-03*T - 8.94616113e-07*pow(T, 2) - 1.49534151e-09*pow(T, 3)	# O2
cp_2[:,2]  = 3.35826314e+01 + 1.46931655e-02*T + 5.99930802e-06*pow(T, 2) + 3.53971581e-09*pow(T, 3)	# SO2
cp_2[:,3]  = 3.01308973e+01 + 8.36012484e-03*T + 1.02417826e-06*pow(T, 2) - 5.99505721e-11*pow(T, 3)	# H2O
cp_2[:,4]  = 3.28975973e+01 + 7.14163891e-03*T - 6.76694887e-07*pow(T, 2) - 1.39803868e-09*pow(T, 3)	# S2
cp_2[:,5]  = 1.95675260e+01 + 4.26884125e-03*T - 3.83532582e-07*pow(T, 2) - 8.17473952e-10*pow(T, 3)	# N2

dH_2[:,0]  = -2.25228053e+04 + 3.67545353e+00*T + 8.20354966e-03*pow(T, 2) + 6.70367713e-06*pow(T, 3)
dH_2[:,1]  = -2.70324897e+03 + 7.85808005e+00*T + 7.68747370e-03*pow(T, 2) + 5.78093859e-06*pow(T, 3)
dH_2[:,2]  = -2.89877255e+05 - 3.50909552e+01*T + 3.11201895e-02*pow(T, 2) + 3.09962707e-05*pow(T, 3)
dH_2[:,3]  = -2.32311516e+05 - 3.82980803e+01*T + 1.57441057e-02*pow(T, 2) + 1.80135000e-05*pow(T, 3)
dH_2[:,4]  = +1.18666906e+05 + 3.23417401e+01*T + 3.50618677e-03*pow(T, 2) - 6.27672665e-07*pow(T, 3)
dH_2[:,5]  = -2.06273070e+03 + 5.03062895e+00*T + 4.99340479e-03*pow(T, 2) + 3.76336698e-06*pow(T, 3)

drH_2[:,0]  = (1.0*dH_2[:,2]  + 1.0*dH_2[:,3] ) - (1.0*dH_2[:,0]  + 3/2*dH_2[:,1])
drH_2[:,1]  = (3/2*dH_2[:,4]  + 2.0*dH_2[:,3] ) - (2.0*dH_2[:,0]  + 1.0*dH_2[:,2])



plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T, cp_1[:,0])
plt.ylim(25, 50)

plt.subplot(2,1,2)
plt.plot(T, cp_2[:,0])
plt.ylim(25, 50)
plt.show()

cp_diff = cp_1[:,0] - cp_2[:,0]

plt.figure(2)
plt.plot(T, cp_diff)
plt.show()


# plt.figure(1)
# plt.subplot(2,1,1)
# plt.plot(T, cp_1[:,0], T, cp_1[:,1], T, cp_1[:,2], T, cp_1[:,3], T, cp_1[:,4], T, cp_1[:,5])
# plt.ylim(25, 85)
# plt.title("Thermal Capacity")
# plt.legend(["H2S", "SO2", "S2", "H2O", "O2", "N2"], loc = "best")
# plt.xlabel("[K]")
# plt.ylabel("[J/molK]")
# plt.subplot(2,1,2)
# plt.plot(T, cp_2[:,0], T, cp_2[:,2], T, cp_2[:,4], T, cp_2[:,3], T, cp_2[:,1], T, cp_2[:,5])
# plt.ylim(25, 85)
# plt.legend(["H2S", "SO2", "S2", "H2O", "O2", "N2"], loc = "best")
# plt.xlabel("[K]")
# plt.ylabel("[J/molK]")
# plt.show()



# plt.title("Concentration Progress")
# plt.xlabel("Time [s]")
# plt.ylabel("Concentration [mol/dm3]")
# plt.legend(["Concentration"], loc = "best")













