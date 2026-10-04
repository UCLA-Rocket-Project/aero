clear; clc;
% to use this, you will need mach number, characteristic length, and altitude
in2m = 0.0254;
feet2m = in2m*12;

%%% Enter here
bodylength = 205 * in2m; % characteristic length of object (i.e. airfoil mean aerodynamic chord or length of rocket)
alt = 2600*feet2m;
mach = 0.8;
r = 1.2; % growth rate for boundary layers. 1.2 is pretty standard
k = 30e-6; % roughness height (5e-5 = 50 microns ~ regular paint)
y_plus = 0.75; % y+ is a nondimensional distance from the wall (relative to freestream)
%%%


[rho, a] = atmo(alt, mach);
vel = mach*a;
function [rho, a] = atmo(alt, mach)
    % Constants
    P0 = 94200;     % air pressure at FAR
    T0 = 298.15;    % 25 degrees celsius, likely for FAR in late spring
    L  = 0.0065;    % Temperature lapse rate (K/m)
    R  = 287.05;    % Gas constant for air (J/kg*K)
    g  = 9.81;
    gamma = 1.4;    % Ratio of specific heats

    % Current Temperature (Standard Lapse Model)
    temp = T0 - L * alt;
    temp = max(temp, 216.65) % Floor at Tropopause temp
    a = sqrt(gamma * R * temp);
    pressure = P0 * (1 - L*alt/T0)^(g/(L*R))
    tot_pressure = pressure*(1+(1.4-1)/2*mach^2)^(1.4/(1.4-1))
    rho = pressure / (R * temp);
end

Re = rho*vel*bodylength / 1.789e-5; 

if(mach<1)
    Re_cutoff = 38.21*(bodylength/(k))^1.053;
else
    Re_cutoff = 44.62*(bodylength/k)^1.053 * mach^1.16;
end
Re_cfd = min(Re, Re_cutoff);
C_f_cfd = 0.455/(log10(Re_cfd)^2.58*(1+0.144*mach^2)^0.65);
first_layer_height = y_plus * 1.789e-5 / rho / sqrt(0.5*vel^2*C_f_cfd)
BL_height = 0.37*bodylength/(Re_cfd^0.2)
N = log(1 - (BL_height * (1 - r) / first_layer_height)) / log(r);
N = ceil(N)
bias = r ^ N
K_plus = sqrt(0.5*vel^2*C_f_cfd) * k / 1.789e-5

turb_I = 0.16*Re^(-1/8)

N2 = log(1 - (0.01016 * (1 - r) / first_layer_height)) / log(r)
bias2 = r ^ N2
% lastLayerHeight = r^(N2+1)
% % N3 = log(1 - (10.98 * (1 - r) / lastLayerHeight)) / log(r)
% % bias3 = r^N3

