%% TIP_THERMO  Stagnation-point heating of the TR1 aluminum nose-cone tip
% Reads an OpenRocket flight-data CSV export and computes:
%   - stagnation-point heat flux q_s'' (Sutton-Graves, cold wall)
%   - integrated heat input per area and total heat into the tip
%   - max heat flux and its time
%   - an ANSYS heat-flux-vs-time table (CSV)
%   - a conservative lumped (no-loss) tip temperature estimate
%
% OpenRocket export (Analyze > Simulation > Export data), any of these columns,
% in any order, any OpenRocket unit:
%   Time, Altitude, Total velocity, Air temperature, Air pressure,
%   Air density, Speed of sound, Mach number
% Minimum needed: Time + (Total velocity or Mach & Speed of sound) + density
% (or pressure & temperature, or altitude).
%
% If "Include field descriptions" is checked on export, the "# Time (s),..."
% header line is used to find columns and units automatically. Otherwise the
% column layout is picked from the `layouts` table below by column count.
%
% Same method as the "thermo" sheet in "TR1 alum tip (fixed).xlsx".

clear; clc; close all;

%% ======================= USER INPUTS =======================
csvFile = 'example ork containing time altitude vel temp density speed of sound mach.csv';

% --- Headerless-export layouts (by number of columns) ---
layouts = {
    8, {'time','alt','vel','temp','pres','rho','a','mach'}
    4, {'time','vel','rho','mach'}
};
% --- Units used when the file has no header (OpenRocket imperial defaults) ---
defUnits = struct('time','s','alt','ft','vel','ft/s','temp','F','pres','mbar', ...
                  'rho','oz/in3','a','ft/s','mach','');

% --- Geometry (in) : Von Karman / LD-Haack tip ---
tipLength_in   = 1.1438;   % tip length
tipBaseR_in    = 0.5;      % tip base radius
noseR_in       = 0.5;      % effective nose radius R_n for Sutton-Graves
tipVolume_in3  = [];       % [] = compute from Haack profile, or enter CAD value
tipArea_in2    = [];       % [] = compute from Haack profile, or enter CAD value

% --- Material: aluminum 7075-T6 ---
rho_al  = 2810;            % density, kg/m^3
cp_al   = 960;             % specific heat, J/kg-K
k_al    = 130;             % thermal conductivity, W/m-K (for ANSYS; not used in lumped calc)
Tmelt_K = 750;             % solidus (start of melting), K  (~477 C)

% --- Heating model ---
kSG      = 1.7415e-4;      % Sutton-Graves constant for air, kg^0.5/m  (q = k*sqrt(rho/Rn)*V^3, W/m^2)
useAirspeed = true;        % true: V = Mach * speed of sound (air-relative); false: OpenRocket total velocity
heatedArea_in2 = [];       % [] = min(2*pi*Rn^2, tip surface area): whole cap at stagnation flux (conservative)

% --- ANSYS output ---
ansysDt   = 0.05;          % resample step for ANSYS table, s ([] = keep raw OpenRocket time steps)
ansysTEnd = [];            % stop the table at this time, s ([] = apogee)

%% ======================= READ FLIGHT DATA =======================
[d, u] = readOpenRocket(csvFile, layouts, defUnits);
fld = @(f) isfield(d,f);

t = d.time(:);
% convert everything that exists to SI
if fld('alt'),  alt = toSI(d.alt,  u.alt);   end
if fld('vel'),  Vtot = toSI(d.vel, u.vel);   end
if fld('temp'), Tinf = toSI(d.temp, u.temp); end
if fld('pres'), pinf = toSI(d.pres, u.pres); end
if fld('rho'),  rho  = toSI(d.rho,  u.rho);  end
if fld('a'),    a    = toSI(d.a,    u.a);    end
if fld('mach'), M    = d.mach(:);            end

gam = 1.4;  Rgas = 287.05;
% fill in anything missing
if ~exist('Tinf','var')
    if exist('a','var'),        Tinf = a.^2/(gam*Rgas);
    elseif exist('alt','var'),  [Tinf,~,~] = isa1976(alt);
    else,                       Tinf = nan(size(t));
    end
end
if ~exist('rho','var')
    if exist('pinf','var') && all(isfinite(Tinf)), rho = pinf./(Rgas*Tinf);
    elseif exist('alt','var'),  [~,~,rho] = isa1976(alt);
    else, error('Need air density, pressure+temperature, or altitude.');
    end
end
if ~exist('a','var') && all(isfinite(Tinf)), a = sqrt(gam*Rgas*Tinf); end
if ~exist('M','var')
    if exist('a','var') && exist('Vtot','var'), M = Vtot./a;
    else, error('Need Mach number, or total velocity plus speed of sound/temperature.'); end
end
if useAirspeed && exist('a','var')
    V = M.*a;                  % air-relative speed
elseif exist('Vtot','var')
    V = abs(Vtot);
else
    error('No velocity available.');
end

%% ======================= GEOMETRY =======================
in = 0.0254;
L  = tipLength_in*in;  R = tipBaseR_in*in;  Rn = noseR_in*in;
[Vh, Ah] = haackVolArea(L, R, 100);
if isempty(tipVolume_in3), Vtip = Vh; else, Vtip = tipVolume_in3*in^3; end
if isempty(tipArea_in2),   Atip = Ah; else, Atip = tipArea_in2*in^2;   end
if isempty(heatedArea_in2), Aheat = min(2*pi*Rn^2, Atip); else, Aheat = heatedArea_in2*in^2; end
mTip = rho_al*Vtip;

%% ======================= HEATING =======================
q   = kSG*sqrt(rho/Rn).*V.^3;                    % W/m^2, stagnation point, cold wall
qInt = cumtrapz(t, q);                           % J/m^2, integrated heat input per area
Qin  = Aheat*qInt;                               % J into tip
T0   = Tinf.*(1 + (gam-1)/2*M.^2);               % stagnation (recovery) temperature, K
Tinit = Tinf(1); if ~isfinite(Tinit), Tinit = 288.15; end
Ttip = Tinit + Qin/(mTip*cp_al);                 % lumped, no re-radiation / conduction losses

[qMax, iMax] = max(q);
iApo = []; if exist('alt','var'), [~, iApo] = max(alt); end

%% ======================= REPORT =======================
fprintf('==== TR1 aluminum tip - stagnation heating ====\n');
fprintf('File: %s  (%d rows, %.2f s)\n', csvFile, numel(t), t(end));
fprintf('Tip: L = %.4f in, R = %.3f in, Rn = %.3f in\n', tipLength_in, tipBaseR_in, noseR_in);
fprintf('     V = %.4f in^3, A = %.4f in^2, heated A = %.4f in^2, m = %.2f g\n', ...
        Vtip/in^3, Atip/in^2, Aheat/in^2, mTip*1e3);
fprintf('Max Mach                 : %.3f  at t = %.2f s\n', max(M), t(find(M==max(M),1)));
fprintf('Max heat flux q_s''''      : %.3f W/cm^2 (%.4g W/m^2)  at t = %.2f s\n', qMax/1e4, qMax, t(iMax));
fprintf('Integrated heat input    : %.2f J/cm^2\n', qInt(end)/1e4);
fprintf('Total heat into tip      : %.1f J\n', Qin(end));
if any(isfinite(T0))
    fprintf('Max stagnation temp T0   : %.1f K (%.1f C)\n', max(T0), max(T0)-273.15);
else
    fprintf('Max stagnation temp T0   : n/a (no temperature/speed of sound/altitude in export)\n');
    fprintf('                           lumped estimate started at 288.15 K\n');
end
fprintf('Lumped tip temp (max)    : %.1f K (%.1f C), rise %.1f K\n', Ttip(end), Ttip(end)-273.15, Ttip(end)-Tinit);
fprintf('Melt margin (solidus)    : %.1f K\n', Tmelt_K - Ttip(end));
if any(isfinite(T0)) && Ttip(end) > max(T0)
    fprintf(['NOTE: lumped estimate exceeds max stagnation temperature; the tip cannot\n' ...
             '      physically get hotter than T0, so the no-loss cold-wall bound is very conservative.\n']);
end

%% ======================= ANSYS TABLE =======================
if isempty(ansysTEnd)
    if ~isempty(iApo), ansysTEnd = t(iApo); else, ansysTEnd = t(find(q >= 1e-3*qMax, 1, 'last')); end  % no altitude: end when q < 0.1% of max
end
if isempty(ansysDt)
    tA = t(t <= ansysTEnd);  qA = q(t <= ansysTEnd);
else
    [tu, iu] = unique(t);
    tA = (0:ansysDt:ansysTEnd)';  qA = interp1(tu, q(iu), tA, 'linear', 0);
end
[~, base] = fileparts(csvFile);
outFile = [base '_ansys_heatflux.csv'];
fid = fopen(outFile, 'w');
fprintf(fid, 'Time [s],Heat Flux [W/m^2]\n');
fprintf(fid, '%.4f,%.6g\n', [tA qA]');
fclose(fid);
fprintf('ANSYS heat-flux table    : %s (%d rows, 0 to %.2f s)\n', outFile, numel(tA), ansysTEnd);

%% ======================= PLOTS =======================
figure('Name','Tip heating','Color','w');
tiledlayout(2,2,'TileSpacing','compact');
nexttile; plot(t, q/1e4,'LineWidth',1.2); hold on; plot(t(iMax), qMax/1e4,'ro');
xlabel('t (s)'); ylabel('q_s'''' (W/cm^2)'); title('Stagnation heat flux'); grid on;
nexttile; plot(t, qInt/1e4,'LineWidth',1.2); xlabel('t (s)'); ylabel('\int q dt (J/cm^2)');
title('Integrated heat input'); grid on;
nexttile; plot(t, Ttip-273.15,'LineWidth',1.2); hold on; plot(t, T0-273.15,'--');
yline(Tmelt_K-273.15,'r:','solidus'); xlabel('t (s)'); ylabel('T (C)');
legend('lumped tip','stagnation T_0','Location','best'); title('Temperature'); grid on;
nexttile; plot(t, M,'LineWidth',1.2); xlabel('t (s)'); ylabel('Mach'); grid on;
if exist('alt','var'), yyaxis right; plot(t, alt/0.3048); ylabel('Altitude (ft)'); end
title('Flight');

%% ======================= FUNCTIONS =======================
function [d, u] = readOpenRocket(file, layouts, defUnits)
% Read numeric data, skipping OpenRocket '#' comment/event lines.
    raw = readmatrix(file, 'CommentStyle', '#', 'FileType', 'text');
    raw = raw(all(isfinite(raw(:,1)),2), :);
    nCol = size(raw, 2);
    % look for a "# Time (s),..." field-description line
    names = {}; units = {};
    txt = splitlines(fileread(file));
    for k = 1:min(numel(txt), 50)
        ln = strtrim(txt{k});
        if startsWith(ln,'#') && contains(lower(ln),'time') && contains(ln,',')
            parts = strsplit(strtrim(ln(2:end)), ',');
            if numel(parts) == nCol
                for p = 1:nCol
                    tok = regexp(parts{p}, '^\s*(.*?)\s*\((.*)\)\s*$', 'tokens', 'once');
                    if isempty(tok), tok = {strtrim(parts{p}), ''}; end
                    names{p} = mapName(tok{1}); units{p} = tok{2}; %#ok<AGROW>
                end
                break
            end
        end
    end
    if isempty(names)
        row = find([layouts{:,1}] == nCol, 1);
        if isempty(row), error('No layout defined for %d columns - add one to "layouts".', nCol); end
        names = layouts{row,2};
        units = cellfun(@(f) defUnits.(f), names, 'UniformOutput', false);
    end
    d = struct(); u = struct();
    for p = 1:nCol
        if isempty(names{p}), continue; end
        d.(names{p}) = raw(:,p); u.(names{p}) = units{p};
    end
    if ~isfield(d,'time'), error('No time column found.'); end
end

function f = mapName(s)
    s = lower(s);
    keys = {'time','time'; 'altitude','alt'; 'total velocity','vel'; 'air temperature','temp';
            'air pressure','pres'; 'air density','rho'; 'speed of sound','a'; 'mach','mach'};
    f = '';
    for k = 1:size(keys,1)
        if contains(s, keys{k,1}), f = keys{k,2}; return; end
    end
end

function x = toSI(x, unit)
    x = x(:);
    switch lower(strrep(strrep(unit,' ',''),'³','3'))
        case {'','s','m','m/s','k','pa','kg/m3'}               % already SI
        case 'ft',      x = x*0.3048;
        case 'in',      x = x*0.0254;
        case 'km',      x = x*1000;
        case 'mi',      x = x*1609.344;
        case 'ft/s',    x = x*0.3048;
        case 'mph',     x = x*0.44704;
        case 'km/h',    x = x/3.6;
        case 'kt',      x = x*0.514444;
        case {'°f','f'},x = (x-32)*5/9 + 273.15;
        case {'°c','c'},x = x + 273.15;
        case 'mbar',    x = x*100;
        case 'bar',     x = x*1e5;
        case 'kpa',     x = x*1e3;
        case 'atm',     x = x*101325;
        case 'psi',     x = x*6894.757;
        case 'inhg',    x = x*3386.389;
        case 'mmhg',    x = x*133.3224;
        case 'oz/in3',  x = x*1729.994;
        case 'lb/in3',  x = x*27679.90;
        case 'lb/ft3',  x = x*16.01846;
        case 'g/cm3',   x = x*1000;
        case 'slug/ft3',x = x*515.3788;
        otherwise, error('Unknown unit "%s" - add it to toSI().', unit);
    end
end

function [V, A] = haackVolArea(L, R, N)
% LD-Haack (Von Karman) profile, theta-spaced points; frustum sums.
    th = linspace(0, pi, N+1);
    x  = L*(1 - cos(th))/2;
    r  = R/sqrt(pi)*sqrt(th - sin(2*th)/2);
    dx = diff(x); dr = diff(r); r1 = r(1:end-1); r2 = r(2:end);
    V = sum(pi/3*dx.*(r1.^2 + r1.*r2 + r2.^2));
    A = sum(pi*(r1 + r2).*hypot(dx, dr));
end

function [T, p, rho] = isa1976(h)
% Standard atmosphere, 0-20 km (geometric ~ geopotential at these heights).
    T = 288.15 - 0.0065*min(h,11000);
    p = 101325*(T/288.15).^5.2559;
    hi = h > 11000;
    p(hi) = 22632*exp(-9.80665*(h(hi)-11000)/(287.05*216.65));
    rho = p./(287.05*T);
end
