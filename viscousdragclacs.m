clear; clc; close all;

%viscous (skin-friction) drag on blunt spike rocket, no fins

%inputs (in)
spk_len = 8;
spk_di = 0.75;
nc_len = 30;
body_len = 153;
bt_len = 15;
aft_di = 5.24;
diam = 8.4;

%wetted area (in^2)
wa_spk_tip = 2*pi*(spk_di/2)^2;        
wa_spk = pi*spk_di*spk_len;
wa_body = pi*diam*body_len;
wa_bt = pi*(diam/2 + aft_di/2)*sqrt(bt_len^2 + (diam/2 - aft_di/2)^2);

r = diam/2;
x = linspace(0,nc_len,1000);
theta = acos(1 - 2*x/nc_len);
y = (r/sqrt(pi))*sqrt(theta - sin(2*theta)/2);   % LD-Haack (Von Karman) nose cone
idx = 2*y >= 0.75;                               % nose cone truncated where spike attaches
x_cut = x(idx);
y_cut = y(idx);
dydx = gradient(y_cut,x_cut);
wa_nc = trapz(x_cut,2*pi*y_cut.*sqrt(1 + dydx.^2));

wa = wa_spk_tip + wa_spk +wa_nc+wa_body+wa_bt;

S_wet = wa * 0.0254^2;                  % in^2 -> m^2
L_rocket = (spk_len + nc_len + body_len + bt_len) * 0.0254;   % m

%surface roughness
k_rough = 1.2e-3;                       % roughness height (in)
k_rough = k_rough * 0.0254;             % in -> m

% ork time v air temp v air density v mach number v reynolds number
data_file = 'flightork time  v air temp v air density v mach number v reynolds number.csv';
col_t   = 1;
col_T   = 2;
col_rho = 3;
col_M   = 4;
col_Re  = 5;


%Gas properties (dry air, calorically perfect)
gam = 1.4;           % ratio of specific heats
R_air = 287.05;      % J/(kg K)
Pr = 0.71;           % Prandtl number
r_rec = Pr^(1/3);    % turbulent recovery factor (~0.89), White, Viscous Fluid Flow

% Sutherland's viscosity law (only the RATIO mu_e/mu* is used)
mu_ref = 1.716e-5;   % Pa s
T_mu_ref = 273.15;   % K
S_suth = 110.4;      % K

% validity lims
Re_min_valid = 5e5;  % Prandtl-Schlichting fit is for turbulent flow, roughly 5e5 < Re < 1e9
Re_max_valid = 1e9;
M_max_valid = 4;     % reference-temperature method checked vs. Van Driest II well past M = 4;

%load open rocket
script_dir = fileparts(mfilename('fullpath'));
data_path = fullfile(script_dir, '..', 'excel misc', data_file);
raw = readmatrix(data_path);

t   = raw(:,col_t);
T   = (raw(:,col_T) - 32)*5/9 + 273.15;         % deg F -> K
rho = raw(:,col_rho) * 0.028349523125/0.0254^3; % oz/in^3 -> kg/m^3 (x 1729.99)
M   = raw(:,col_M);
Re  = raw(:,col_Re);                            % OpenRocket Re, based on total rocket length

%find apogee: last point before Mach starts rising again after peak Mach
[~, iMmax] = max(M);
i_apo = iMmax + find(diff(M(iMmax:end)) > 0, 1) - 1;
t = t(1:i_apo); rho = rho(1:i_apo); M = M(1:i_apo); Re = Re(1:i_apo); T = T(1:i_apo);

% Airspeed from Mach (includes wind, consistent with OpenRocket's Re).
a = sqrt(gam*R_air*T);               % speed of sound (m/s)
V = M .* a;                          % airspeed (m/s)

% SKIN RFRITION COEFF: text is from chat
% Incompressible turbulent flat plate (Prandtl-Schlichting):
%   Cf_inc = 0.455 / (log10(Re_L))^2.58
% Source: Schlichting, Boundary-Layer Theory; White, Viscous Fluid Flow.
% Average Cf over a plate of length L, fully turbulent from x = 0.
%
% Compressibility: Eckert reference-temperature method
%   (Eckert, J. Aeronautical Sciences 22(8), 1955; White, Viscous Fluid Flow, Ch. 7)
%   The incompressible formula is evaluated with gas properties at a
%   reference temperature T* inside the boundary layer:
%     T_aw = T_e*(1 + r*(gam-1)/2*M^2)                 adiabatic-wall temperature
%     T*   = T_e + 0.5*(T_w - T_e) + 0.22*(T_aw - T_e)
%     Re*  = Re_e * (rho*/rho_e)*(mu_e/mu*) = Re_e * (T_e/T*)*(mu_e/mu*)
%     Cf*  = 0.455/(log10 Re*)^2.58
%     Cf   = Cf* * (rho*/rho_e) = Cf* * (T_e/T*)       (referenced to freestream q)
%   (pressure is constant across the boundary layer, so rho ~ 1/T)
% Assumptions: flat plate, zero pressure gradient, perfect gas, Sutherland
%   viscosity, the rocket body treated as a flat plate of length L (body
%   transverse curvature neglected; valid since boundary layer << body radius).
% Range: subsonic through supersonic; agrees with Van Driest II within a few
%   percent for adiabatic walls up to at least M ~ 4-5 (Hopkins & Inouye 1971).
%   At M -> 0 it reduces exactly to the incompressible formula.
%
% Cross-check (not used for the result): Raymer's adiabatic-wall fit
%   Cf = Cf_inc / (1 + 0.144*M^2)^0.65   (Raymer, Aircraft Design, eq. 12.27)
%
% Transonic note: these flat-plate corrections are smooth through M = 1 (skin
% friction depends on M only through the temperature rise in the boundary
% layer). The real error near M = 1 comes from shocks on the nose/boattail
% shoulders and shock / boundary-layer interaction, which can locally
% thicken or separate the boundary layer, so results for 0.8 <= M <= 1.2
% are less reliable.

% Adiabatic wall: Tw = T_aw (insulated, thermally equilibrated skin), the
% standard assumption behind published supersonic skin-friction correlations.
% Real skin lags behind T_aw and stays cooler, which gives a slightly HIGHER Cf.
T_aw = T .* (1 + r_rec*(gam - 1)/2 .* M.^2);
Tw = T_aw;

T_star  = T + 0.5*(Tw - T) + 0.22*(T_aw - T);
mu_e    = mu_ref*(T/T_mu_ref).^1.5*(T_mu_ref + S_suth)./(T + S_suth);   % freestream viscosity (Sutherland)
mu_star =mu_ref*(T_star/T_mu_ref).^1.5*(T_mu_ref + S_suth)./(T_star + S_suth);
Re_star = Re .* (T./T_star) .* (mu_e./mu_star);

% currently asumming no laminar to turblenet transitions but theretically
% flow is laminar on nc
Cf_star = 0.455 ./ (log10(Re_star)).^2.58;
Cf_smooth = Cf_star .* (T./T_star);

% Comparison values (not used for drag)
Cf_inc    = 0.455 ./ (log10(Re)).^2.58;
Cf_raymer = Cf_inc ./ (1 + 0.144*M.^2).^0.65;

% Surface roughness: cutoff Reynolds number (Raymer, Aircraft Design, eqs.
% 12.28-12.29, from USAF DATCOM; same approach as OpenRocket). Above Re_cut
% the roughness elements stick out of the viscous sublayer and Cf stops
% falling with Re, so Cf is evaluated at Re_cut instead of the flight Re.
%   subsonic:              Re_cut = 38.21*(L/k)^1.053
%   transonic/supersonic:  Re_cut = 44.62*(L/k)^1.053*M^1.16
% The larger of the two is used, which switches smoothly at M ~ 0.875 where
% they are equal. The rough Cf uses Raymer's compressibility factor, which
% matches the Eckert smooth result within ~2% over this flight.
% Cf = max(smooth, rough): the rough value only governs when Re > Re_cut.
Re_cut = max(38.21*(L_rocket/k_rough)^1.053, 44.62*(L_rocket/k_rough)^1.053 * M.^1.16);
Cf_rough = 0.455 ./ (log10(Re_cut)).^2.58 ./ (1 + 0.144*M.^2).^0.65;
Cf = max(Cf_smooth, Cf_rough);
rough_governs = Cf_rough > Cf_smooth;

% Only use rows with positive data inside the correlation's valid range (no
% extrapolation). Missing/NaN values fail these comparisons automatically.
use = rho > 0 & M > 0 & M <= M_max_valid ...
    & Re_star >= Re_min_valid & Re_star <= Re_max_valid & Re_cut >= Re_min_valid;
if ~all(use)
    warning('%d points with bad data or outside the valid range were excluded.', nnz(~use));
end
Cf(~use) = NaN;

% dynamic pressure eq
q = 0.5 * rho .* V.^2;           % Pa
%skin drag eq
D_f = q .* Cf * S_wet;           % N

% skin-friction drag coefficient, referenced to body cross-section area
% (same reference area OpenRocket uses): Cd_f = D_f/(q*A_ref) = Cf*S_wet/A_ref
A_ref = pi*(diam*0.0254)^2/4;    % m^2
Cd_f = Cf * S_wet / A_ref;

% Drag impulse
I_total = trapz(t(use), D_f(use));

%Results
fprintf('Boundary layer: fully turbulent, wall: adiabatic, roughness: %.1e in\n', k_rough/0.0254);
fprintf('Skin-friction drag impulse, launch to apogee (0 - %.1f s): %.1f N s\n', t(end), I_total);
fprintf('  (smooth surface would give %.1f N s)\n', trapz(t(use), q(use).*Cf_smooth(use)*S_wet));

%plot
%drag v time
figure('Name','Skin-friction drag');
h1 = plot(t, D_f);
xlabel('Time (s)'); ylabel('Skin-friction drag D_f (N)');
title(sprintf('Skin-friction drag to apogee (fully turbulent, adiabatic wall, k = %.1e in), I = %.0f N s', k_rough/0.0254, I_total));

% Cd v Mach
[~, iMmax] = max(M);
figure('Name','Cd vs Mach');
plot(M(1:iMmax), Cd_f(1:iMmax), M(iMmax:end), Cd_f(iMmax:end));
xlabel('Mach'); ylabel('Skin-friction drag coefficient C_{D,f}');
title(sprintf('Skin-friction C_D vs Mach (A_{ref} = %.4f m^2, %.1f in dia)', A_ref, diam));
legend('Accelerating (to peak Mach)', 'Decelerating (to apogee)');
