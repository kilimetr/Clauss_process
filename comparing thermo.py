# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: Comparing thermo between two models

import numpy as np
import matplotlib.pyplot as plt

from regres import results, dHstd



T = np.linspace(298.15, 1500, 800)



# PYTHON
cp_py    = np.empty([len(T), 6])
dH_py    = np.empty([len(T), 6])
drH_1_py = np.empty([len(T), 1])
drH_2_py = np.empty([len(T), 1])

cp_py[:, 0] = results["H2S"][0][3] + results["H2S"][0][2]*T + results["H2S"][0][1]*T*T + results["H2S"][0][0]*T*T*T   # H2S
cp_py[:, 1] = results["SO2"][0][3] + results["SO2"][0][2]*T + results["SO2"][0][1]*T*T + results["SO2"][0][0]*T*T*T   # SO2
cp_py[:, 2] = results["S2"][0][3]  + results["S2"][0][2] *T + results["S2"][0][1] *T*T + results["S2"][0][0] *T*T*T   # S2
cp_py[:, 3] = results["H2O"][0][3] + results["H2O"][0][2]*T + results["H2O"][0][1]*T*T + results["H2O"][0][0]*T*T*T   # H2O
cp_py[:, 4] = results["O2"][0][3]  + results["O2"][0][2] *T + results["O2"][0][1] *T*T + results["O2"][0][0] *T*T*T   # O2
cp_py[:, 5] = results["N2"][0][3]  + results["N2"][0][2] *T + results["N2"][0][1] *T*T + results["N2"][0][0] *T*T*T   # N2

dH_py[:, 0] = results["H2S"][1][3] + results["H2S"][1][2]*T + results["H2S"][1][1]*T*T + results["H2S"][1][0]*T*T*T   # H2S
dH_py[:, 1] = results["SO2"][1][3] + results["SO2"][1][2]*T + results["SO2"][1][1]*T*T + results["SO2"][1][0]*T*T*T   # SO2
dH_py[:, 2] = results["S2"][1][3]  + results["S2"][1][2] *T + results["S2"][1][1] *T*T + results["S2"][1][0] *T*T*T   # S2
dH_py[:, 3] = results["H2O"][1][3] + results["H2O"][1][2]*T + results["H2O"][1][1]*T*T + results["H2O"][1][0]*T*T*T   # H2O
dH_py[:, 4] = results["O2"][1][3]  + results["O2"][1][2] *T + results["O2"][1][1] *T*T + results["O2"][1][0] *T*T*T   # O2
dH_py[:, 5] = results["N2"][1][3]  + results["N2"][1][2] *T + results["N2"][1][1] *T*T + results["N2"][1][0] *T*T*T   # N2

dH_py[:, 0] = dH_py[:, 0] + dHstd["H2S"]
dH_py[:, 1] = dH_py[:, 1] + dHstd["SO2"]
dH_py[:, 2] = dH_py[:, 2] + dHstd["S2"]
dH_py[:, 3] = dH_py[:, 3] + dHstd["H2O"]
dH_py[:, 4] = dH_py[:, 4] + dHstd["O2"]
dH_py[:, 5] = dH_py[:, 5] + dHstd["N2"]

dH_py[:, 0] = dH_py[:, 0] * 1000
dH_py[:, 1] = dH_py[:, 1] * 1000
dH_py[:, 2] = dH_py[:, 2] * 1000
dH_py[:, 3] = dH_py[:, 3] * 1000
dH_py[:, 4] = dH_py[:, 4] * 1000
dH_py[:, 5] = dH_py[:, 5] * 1000

print(dH_py.shape)

drH_py_1 =     dH_py[:, 1] +   dH_py[:, 3] - (  dH_py[:, 0] + 3/2*dH_py[:,4])
drH_py_2 = 3/2*dH_py[:, 2] + 2*dH_py[:, 3] - (2*dH_py[:, 0] +     dH_py[:,1])

drHstd_py_2 = 3/2*dHstd["S2"] + 2*dHstd["H2O"] - (2*dHstd["H2S"] + dHstd["SO2"])
drHstd_py_2 = drHstd_py_2 * 1000

print(drHstd_py_2)

dcp_1 =     cp_py[:, 1] +   cp_py[:, 3] - (  cp_py[:, 0] + 3/2*cp_py[:,4])
dcp_2 = 3/2*cp_py[:, 2] + 2*cp_py[:, 3] - (2*cp_py[:, 0] +     cp_py[:,1])

# drH_py_11 = 
drH_py_22 = drHstd_py_2 + dcp_2


# MATLAB
cp_m  = np.empty([len(T), 6])
drH_m = np.empty([len(T), 2])

#	FIRST MODEL'S DATA
#			   H2S + 3/2 O2  ->    SO2 +   H2O
#			 2 H2S +    SO2 <=> 3/2 S2 + 2 H2O
#			   Positions in vector: H2S, SO2, S2, H2O, O2, N2


cp_m[:, 0] = 0.001436*T + 0.00002432*pow(T, 2) - 1.176e-8*pow(T, 3) + 31.94		# H2S
cp_m[:, 1] = 0.000845*T - 0.0000431 /pow(T, 2)                      + 36.162	# SO2 - bullshit
cp_m[:, 2] =                                                        + 35.0		# S2
cp_m[:, 3] = 0.0113  *T                                             + 30.12		# H2O
cp_m[:, 4] = 0.007171*T - 0.00008556/pow(T, 2)                      + 47.7		# O2 - bullshit
cp_m[:, 5] = 0.00418 *T                                             + 27.82		# N2 - bullshit

drH_m[:,0] = -(2.31198e-35*(2.19016e40*T + 1.60936e36*pow(T, 2) + 1.02726e30*pow(T, 3) + 3.50637e29*pow(T, 4) - 1.27164e26*pow(T, 5) + 3.68688e30))/T
drH_m[:,1] =  (4.62396e-35*(9.48365e38*T + 2.74613e35*pow(T, 2) + 2.04186e32*pow(T, 3) - 3.50637e29*pow(T, 4) + 1.27164e26*pow(T, 5) - 9.32101e29))/T



plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T, drH_m[:, 0], T, drH_py_1)
plt.legend(["matlab", "python"])

plt.subplot(2,1,2)
plt.plot(T, drH_m[:, 1], T, drH_py_2, T, drH_py_22)
plt.legend(["matlab", "python"])
plt.show()

cp_diff = cp_m[:,0] - cp_py[:,0]
cp_diff = sum(cp_diff)


