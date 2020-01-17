# author Dominik Capkovic 
# contact: domcapkovic@gmail.com; https://www.linkedin.com/in/dominik-čapkovič-b0ab8575/
# GitHub: https://github.com/kilimetr


import numpy as np
import scipy as sc
import os

os.chdir("/Users/kilimetr/Desktop/python/Clauss_process")

# BURNER

p     = 101325
n0H2S = 3
n0SO2 = 0
n0S2  = 0
n0H2O = 0
n0O2  = 3/2
n0N2  = 79/21*n0O2
Tin   = 200+273.15

pars = [p]

yini = [n0H2S n0SO2 n0S2 n0H2O n0O2 n0N2 Tin]

[VV,yy] = ode15s(@(V,y) claus_burner(V,pars,y),[0 1], yini)

figure(1)
subplot(2,1,1)
plot(VV,yy(:,1),VV,yy(:,2),VV,yy(:,3),VV,yy(:,4),VV,yy(:,5),VV,yy(:,6))
ylim([0 ceil(n0N2)])
ylabel('Mol flow [mole/s]'); xlabel('Volume [m3]'); title('Mole Changeover in Volume of Burner');
legend('H2S','SO2','S2','H2O','O2','N2');

subplot(2,1,2);
yy(:,end) = yy(:,end) - 273.15; plot(VV,yy(:,end)); ylabel('Temperature [°C]'); xlabel('Volume [m3]');
title('Temperature Changeover in Volume of Burner'); legend('T [°C]'); box on; grid on;

disp(yy(end,5)); % oxygen leftover



# CATALYTIC REACTOR

n0H2S = yy(end,1); % [mol/s]
n0SO2 = yy(end,2);
n0H2O = yy(end,4);
p     = 101325; % [Pa]
Tin   = 400+273.15;
n0S2  = 0;
n0S6  = 0;
n0S8  = 0;
n0N2  = yy(end,6);
n0O2  = yy(end,5);

pars = [p Tin];

yini = [n0H2S n0SO2 n0S2 n0S6 n0S8 n0H2O n0N2 n0O2 Tin];

[VVV,yyy] = ode15s(@(V,y) claus_cat(V,pars,y),[0 0.000000002], yini);

yyy(:,end) = yyy(:,end)-273.15;

figure(2); 
subplot(2,2,1); 
plot(VVV,yyy(:,1),VVV,yyy(:,2),VVV,yyy(:,3),VVV,yyy(:,4),VVV,yyy(:,5),VVV,yyy(:,6),VVV,yyy(:,7),...
     VVV,yyy(:,8)); xlim([0 0.000000000003]);
grid on; box on; xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); 
title('Mole Changeover in Volume of Catalytic Reactor'); 
legend('nH2S','nSO2','nS2', 'nS6','nS8','nH2O','nN2','nO2');

subplot(2,2,2); plot(VVV,yyy(:,end)); xlabel('Volume [m3]'); ylabel('Temperature [°C]');
title('Temperature Changeover in Volume of Catalytic Reactor'); box on; grid on;

subplot(2,2,3); plot(VVV,yyy(:,1),VVV,yyy(:,2),VVV,yyy(:,3),VVV,yyy(:,4),VVV,yyy(:,5),VVV,yyy(:,6)); 
legend('H2S','SO2','S2','S6','S8','H2O'); box on; grid on; xlabel('Volume [m3]');
xlim([0 0.000000000003]); ylabel('Mole flow [mole/s]'); title('Components Changeover in Reactions');

subplot(2,2,4); plot(VVV,yyy(:,3),VVV,yyy(:,4),VVV,yyy(:,5)); legend('S2','S6','S8'); box on; grid on;
xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); title('Various Sulfur Components');

figure(3);
subplot(2,2,1);
plot(VVV,yyy(:,1),VVV,yyy(:,8),VVV,yyy(:,2),VVV,yyy(:,6)); title('H2S + 3/2O2 <=> SO2 + H2O'); box on; 
grid on; xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); legend('H2S','O2','SO2','H2O');
xlim([0 0.000000000003]);

subplot(2,2,2);
plot(VVV,yyy(:,1),VVV,yyy(:,2),VVV,yyy(:,3),VVV,yyy(:,6)); title('2H2S + SO2 <=> 3/2S2 + 2H2O'); box on; 
grid on; xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); legend('H2S','SO2','S2','H2O');
xlim([0 0.000000000003]);

subplot(2,2,3);
plot(VVV,yyy(:,1),VVV,yyy(:,2),VVV,yyy(:,4),VVV,yyy(:,6)); title('2H2S + SO2 <=> 3/6S6 + 2H2O'); box on; 
grid on; xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); legend('H2S','SO2','S6','H2O');
xlim([0 0.000000000003]);

subplot(2,2,4);
plot(VVV,yyy(:,1),VVV,yyy(:,2),VVV,yyy(:,5),VVV,yyy(:,6)); title('2H2S + SO2 <=> 3/8S8 + 2H2O'); box on; 
grid on; xlabel('Volume [m3]'); ylabel('Mole flow [mole/s]'); legend('H2S','SO2','S8','H2O');
xlim([0 0.000000000003]);

XSO2 = (n0SO2 - yy(end,2)) / n0SO2;


return;
%%% thermostuff
Tin     = 298; % K
Tmax    = 400+273.15;
Tvar    = zeros(1,200);
Tvar(1) = Tin;
dT      = (Tmax-Tin) / length(Tvar);

dHH2O = -70000;

for i=2:length(Tvar)
    Tvar(i) = Tvar(i-1)+dT;
end


% kinetic dependence on temp
k0S   = 6.07*10^(-5);
ES    = 3.08*10^(+4);
K0H2O = 0.34*10^(-3);
R     = 8.314;

for i=1:length(Tvar)
    kS(i)   = k0S * exp(-ES/(R*Tvar(i)));               % [molH2S/kg/s]
    KH2O(i) = K0H2O * exp(-dHH2O/(R*Tvar(i)));          % [1/Pa]
    KE(i)   = 9.502*10^(-7) * exp(1.11*10^(4)/Tvar(i)); % [1/kPa^0.5]
    KK(i)    = KE(i)^(0.5)/1000;                         % [1/Pa^0.25]
end

figure(4);
subplot(2,2,1); plot(Tvar,kS);   ylabel('kS [mol "i"/kg/s/Pa^1.5]'); xlabel('T [K]');
subplot(2,2,2); plot(Tvar,KH2O); ylabel('KH2O [1/Pa]');              xlabel('T [K]');
subplot(2,2,3); plot(Tvar,KE);   ylabel('KE [1/kPa^0.5]');           xlabel('T [K]');
subplot(2,2,4); plot(Tvar,KK);    ylabel('K [1/Pa^0.25]');            xlabel('T [K]');

% equilibrium
dHH2S = zeros(length(Tvar),1);
dHSO2 = zeros(length(Tvar),1);
dHS6  = zeros(length(Tvar),1);
dHH2O = zeros(length(Tvar),1);

K = zeros(length(Tvar),1);

dHH2S(1) =  -21000; % [J/mol]
dHSO2(1) = -296830;
dHS6(1)  =  101922;
dHH2O(1) = -241818;
dHS2     =  128600;
dHS8     =  100416;

dSH2S = 205.69; % [J/K/mol]
dSSO2 = 248.11;
dSS6  = 228.16;
dSH2O = 188.72;
dSS2  = 228.16;
dSS8  =  32.04;

% dGH2S = -33400; dGSO2 = -300190; dGS6 = 53699; dGH2O = -228589; dGS2 = 79687; dGS8 = 48578;

for i=1:length(Tvar)
    T = Tvar(i);
    
    cpH2S(i) = 31.94  + (0.1436e-2)*T + (0.2432e-4)*T^2 + (-1.176e-8)*T^3; % [J/K/mol]
    cpSO2(i) = 36.162 +  0.845e-3  *T -    4.310e-5/T^2;
    cpS6(i)  = 127;
    cpH2O(i) = 30.120 +  11.300e-3 *T;
    cpS2(i)  = 35     +           0*T;
    cpS8(i)  = 175;
    
    drcp(i) = 0.5*cpS6(i) + 2*cpH2O(i) - 2*cpH2S(i) - 1*cpSO2(i); % reaction heat capacity
end

for i=2:length(Tvar)
    T = Tvar(i);
    
    dHH2S(i) = dHH2S(1) + cpH2S(i)*(T-293.15);
    dHSO2(i) = dHSO2(1) + cpSO2(i)*(T-293.15);
    dHS6(i)  = dHS6(1)  + cpS6(i) *(T-293.15);
    dHH2O(i) = dHH2O(1) + cpH2O(i)*(T-293.15);
end

for i=1:length(Tvar)
    drH(i) = 0.5*dHS6(i) + 2*dHH2O(i) - 2*dHH2S(i) - 1*dHSO2(i); % reaction enthalpy
    drS    = 0.5*dSS6    + 2*dSH2O    - 2*dSH2S    - 1*dSSO2;    % reaction enthropy
    
    drG(i) = drH(i) - Tvar(i)*drS;
    
    K(i) = exp(-drG(i)/(R*Tvar(i)));
end

figure(5);
subplot(2,3,1); plot(Tvar,cpH2S); xlabel('Temp [K]'); ylabel('cPH2S [J/K/mol]');
subplot(2,3,2); plot(Tvar,cpSO2); xlabel('Temp [K]'); ylabel('cpSO2 [J/K/mol]');
subplot(2,3,3); plot(Tvar,cpS6);  xlabel('Temp [K]'); ylabel('cpS6 [J/K/mol]');
subplot(2,3,4); plot(Tvar,cpH2O); xlabel('Temp [K]'); ylabel('cpH2O [J/K/mol]');
subplot(2,3,5); plot(Tvar,cpS2);  xlabel('Temp [K]'); ylabel('cpS2 [J/K/mol]');
subplot(2,3,6); plot(Tvar,cpS8);  xlabel('Temp [K]'); ylabel('cpS8 [J/K/mol]');

figure(6);
subplot(2,2,1); plot(Tvar,dHH2S); xlabel('Temp [K]'); ylabel('Enthalpy of H2S [J/mol]');
subplot(2,2,2); plot(Tvar,dHSO2); xlabel('Temp [K]'); ylabel('Enthalpy of SO2 [J/mol]');
subplot(2,2,3); plot(Tvar,dHS6);  xlabel('Temp [K]'); ylabel('Enthalpy of S6 [J/mol]');
subplot(2,2,4); plot(Tvar,dHH2O); xlabel('Temp [K]'); ylabel('Enthalpy of H2O [J/mol]');


figure(7);
subplot(2,2,1); plot(Tvar,drcp); xlabel('Temp [K]'); ylabel('Reaction Heat Capacity [J/K/mol]');
subplot(2,2,2); plot(Tvar,drH);  xlabel('Temp [K]'); ylabel('Reaction Enthalpy [J/mol]');
subplot(2,2,3); plot(Tvar,drG);  xlabel('Temp [K]'); ylabel('Reaction Gibbs Energy [J/mol]');

figure(8);
plot(Tvar,K); xlabel('Temp [K]'); ylabel('Equilibrium constant [-]');

