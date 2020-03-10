# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr
# Description: Fit raw data of C2H6 & C3H8
# cp    heat capacity    [J/mol/K]
# dH    enthalpy		 [J/mol]

import numpy as np 
import matplotlib.pyplot as plt 


T_C2H6_raw  = [100, 200, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600]
cp_C2H6_raw = [35.70, 42.30, 52.49, 52.71, 65.46, 77.94, 89.19, 99.14, 107.94, 115.71, 122.55, 128.55, 133.80, 138.39, 142.40, 145.90, 148.98, 151.67, 154.04, 156.14,
			  158.00, 159.65, 161.12, 162.43, 163.61, 164.67, 165.63]
dH_C2H6_std = -83800

T_C3H8_raw  = [50, 100, 150, 200, 273.15, 298.15, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]
cp_C3H8_raw = [34.06, 41.30, 48.79, 56.07, 68.74, 73.60, 73.93, 94.01, 112.59, 128.70, 142.67, 154.77, 165.35, 174.60, 182.67, 189.74, 195.85, 201.21 ,205.89]
dH_C3H8_std = -104700

cp_C2H6_coef = np.polyfit(T_C2H6_raw, cp_C2H6_raw, 3)
cp_C3H8_coef = np.polyfit(T_C3H8_raw, cp_C3H8_raw, 3)

T = np.linspace(100, 1500, 10)

cp_C2H6_fit = cp_C2H6_coef[3] + cp_C2H6_coef[2]*T + cp_C2H6_coef[1]*pow(T, 2) + cp_C2H6_coef[0]*pow(T, 3)
cp_C3H8_fit = cp_C3H8_coef[3] + cp_C3H8_coef[2]*T + cp_C3H8_coef[1]*pow(T, 2) + cp_C3H8_coef[0]*pow(T, 3)

dH_C2H6_calc = dH_C2H6_std + cp_C2H6_coef[3]*T + cp_C2H6_coef[2]*pow(T, 2)/2 + cp_C2H6_coef[1]*pow(T, 3)/3 + cp_C2H6_coef[0]*pow(T, 4)/4
dH_C3H8_calc = dH_C3H8_std + cp_C3H8_coef[3]*T + cp_C3H8_coef[2]*pow(T, 2)/2 + cp_C3H8_coef[1]*pow(T, 3)/3 + cp_C3H8_coef[0]*pow(T, 4)/4

dH_C2H6_coef = np.polyfit(T, dH_C2H6_calc, 3)
dH_C3H8_coef = np.polyfit(T, dH_C3H8_calc, 3)

dH_C2H6_fit = dH_C2H6_coef[3] + dH_C2H6_coef[2]*T + dH_C2H6_coef[1]*pow(T, 2) + dH_C2H6_coef[0]*pow(T, 3)
dH_C3H8_fit = dH_C3H8_coef[3] + dH_C3H8_coef[2]*T + dH_C3H8_coef[1]*pow(T, 2) + dH_C3H8_coef[0]*pow(T, 3)

print(cp_C2H6_coef)
print(cp_C3H8_coef)
print(dH_C2H6_coef)
print(dH_C3H8_coef)

plt.figure(1)
plt.subplot(2, 2, 1)
plt.plot(T_C2H6_raw, cp_C2H6_raw)

plt.subplot(2, 2, 2)
plt.plot(T_C3H8_raw, cp_C3H8_raw)

plt.subplot(2, 2, 3)
plt.plot(T, cp_C2H6_fit)

plt.subplot(2, 2, 4)
plt.plot(T, cp_C3H8_fit)
plt.show()


plt.figure(2)
plt.subplot(2, 2, 1)
plt.plot(T, dH_C2H6_fit)

plt.subplot(2, 2, 2)
plt.plot(T, dH_C3H8_fit)

plt.subplot(2, 2, 3)
plt.plot(T, dH_C2H6_fit)

plt.subplot(2, 2, 4)
plt.plot(T, dH_C3H8_fit)
plt.show()
