% RADA (Rapid Apogee Design Analysis) - by Logan Trickey
% version 3.0
% references: Mandell, Topics in Advanced Model Rocketry
% Hoerner, Fluid-Dynamic Drag

% ~~~~~~INSTRUCTIONS~~~~~~: 
% make sure you have the thrust curve file imported as your
% file path won't match mine. Edit the variable 'mode' to whatever
% parameter you wish to optimize for apogee (parameter names are listed
% there as well). Also change the upper and lower bounds for the parameter
% you are optimizing (avoid unrealistic values like 25 inches for
% semi_span, and select an appropriate step_count value (more step count
% gives a smoother graph at the cost of computing power; I recommend 20-30)

% Version 3 patch notes:
% fixed forebody drag formula, now corroborating RasAero predictions

% ~~~~~~FEATURES~~~~~~:
% 3.A: Now includes option to output the entire damping ratio profile for a
% single flight. To enable, simply set plotDampingProfile from 0 --> 1 and
% set a target_semi_span (when operating with mode = semi_span, the profile
% will plot the damping for semi_span = target_semi_span only)

% 3.B: Now includes option to obtain drag + mass impulse from RasAero for
% validation. To enable, simply set getRasDragImpulse to 1 and upload a
% csv file from rasaero, exporting Mach vs. Cd. Total Impulse will be
% outputted in the command window.

% ACKNOWLEDGEMENTS
% Keep in mind this only accounts for drag and weight
% forces along with thrust, so just because you get a high apogee doesn't
% mean the rocket's stability will be sound. Ultimately, you still need to
% find a balance between static and dynamic stability, damping ratio, and
% the potential for failure modes such as buckling, flutter, etc. The
% reason we made this when OpenRocket or RasAero have a similar function is
% because they calculate mass incorrectly and this also allows you to run
% as many tests as you want at once
% Another limitation is that boattail optimization doesn't really work 
% (unable to find equations that account for flow seperation)

% TERMS AND CONDITIONS
% By utilizing this script you confirm that you have agreed to the user
% license agreement and that you comply with its terms. By utilizing this
% script you agree to sell your soul to the aero gods and submit to
% death in response to any violations of Bernoulli's laws.

clear; clc;
%%
% unit conversions
global in_to_m 
in_to_m = 0.0254;
global lbs_to_kg 
lbs_to_kg = 0.45359237;
g = 9.80665;

plotApogee = 0;
plotDamping = 1;
plotCompMass = 0;
target_semi_span = 7 * in_to_m;         % Semi-span in inches for single test
plotDampingProfile = 1;       % Plot damping vs time for single flight
fintype = 2; % 1 = foam core, 2 = alum


% INPUT HERE YOUR ENGINE PARAMETERS
thrust_curve = readmatrix("C:\Users\ltric\Downloads\Ares 26-27\Thrust Curves\thrustcurvetest2.csv", NumHeaderLines=1);
thrust_time = thrust_curve(:, 1);
thrust_force = thrust_curve(:, 2);
burn_time = 14.9;
prop_mass = 41.7 * lbs_to_kg;
len_tank = 130 * in_to_m;

% ALSO ENTER
rocket_mass_no_aero = 242*lbs_to_kg;
ref_diam = 8.2 * in_to_m;
len_bodytubes = 163*in_to_m;
CG_no_fins = 48.28 * in_to_m;
height_shoulder = 4*in_to_m;



% ENTER YOUR LAYER COUNTS FOR FIBERGLASS AND CF
fin_layers = 6;
nc_layers = 5;
nc_shoulder_layers = 5;
bt_layers = 5;
fin_thickness = 0.25 * in_to_m; % only need to set if you are using aluminum fins, otherwise, check fin_layers


rho_fiberglass = 2000; % Density of 9oz fiberglass [kg/m³] based on an estimate of mass/volume of our test rocket nose
rho_cf = 1560;
rho_foam = 75;
rho_alum = 2700;
thick_cf = 0.0083 * in_to_m;
thick_fg = 0.0003048;
        

%%


getRasAeroDragImpulse = 0; % enable this to obtain nose cone total impulse (mass + drag impulse) using RasAero Mach vs Cd plot.
if getRasAeroDragImpulse
    file = "C:\Users\ltric\Downloads\Ares 26-27\CD TestJT.CSV";
    rasdata = readmatrix(file);
    rasmach = rasdata(:,1);
    rascd = rasdata(:,3);
else
    rasmach = 0; rascd = 0;
end


% enter parameters for optimization

% options for mode: 'len_nose', 'bt_r2', 'len_bt',
% 'semi_span', 'c' (root chord), 'tip' (tip chord)
% 'sweep' (sweep length)

mode = 'semi_span';
lower_bound = 7 * in_to_m;
upper_bound = 55 * in_to_m;
step_count = 1;
step_size = (upper_bound - lower_bound) / step_count;
apogees = zeros(size(step_count));
masses = zeros(size(step_count));
mass_totals = zeros(size(step_count));

max_Q = 0;

if strcmp(mode, 'len_nose')
    params.len_nose = lower_bound;
else
    params.len_nose = 40 * in_to_m;
end

if strcmp(mode, 'len_bt')
    params.len_bt = lower_bound;
else
    params.len_bt = 11 * in_to_m;
end

if strcmp(mode, 'bt_r2')
    params.bt_r2 = lower_bound;
else
    params.bt_r2 = 5.24/2*in_to_m;
end

if strcmp(mode, 'semi_span')
    params.semi_span = lower_bound;
else
    params.semi_span = 7 * in_to_m;
end

if strcmp(mode, 'c')
    params.c = lower_bound;
else
    params.c = 18 * in_to_m;
end

if strcmp(mode, 'tip')
    params.tip = lower_bound;
else
    params.tip = 4 * in_to_m;
end

if strcmp(mode, 'sweep')
    params.sweep = lower_bound;
else
    params.sweep = params.c-params.tip;
end


for i = 1: step_count

    % Parametrizations of mass 
    % ------------- NOSECONE -----------------
    if strcmp(mode, 'len_nose') || i == 1
        len_nose = params.len_nose;
        R = ref_diam / 2;
        A_ref_nose = pi * (ref_diam/2)^2;     % Reference area [m^2]                
        wall_thickness = nc_layers*thick_fg;  % Wall thickness of the nose cone [m]
        g = 9.81;                  % graviational acceleration[m/s^2]
        n_slices = 10000;
        x = linspace(0, len_nose, n_slices);
        % Calculate radius at each position using Von Karman equation
        theta = acos(1 - (2*x)/len_nose);
        y = (R/sqrt(pi)) * sqrt(theta - sin(2*theta)/2);
        r_outer = y;

        dydx = gradient(y, x);
        S_nose = 2 * pi * trapz(x, r_outer.*sqrt(1+dydx.^2));
        r_inner = max(0, y - wall_thickness); % Don't go negative near tip
        V_nosecone = trapz(x, pi * (r_outer.^2 - r_inner.^2));
        V_wave_haack = trapz(x, pi * r_outer.^2);
        %calculate mass of shoulder
        shoulder_wall_thickness = nc_shoulder_layers * 0.000254;
        shoulder_r_outer = ref_diam/2 - 0.003048;
        shoulder_r_inner = shoulder_r_outer - shoulder_wall_thickness;
        shoulder_cross_section = pi * (shoulder_r_outer^2 - shoulder_r_inner^2);
        
        V_shoulder = height_shoulder * shoulder_cross_section;
        %V_shoulder + V_nosecone
        mass_nosecone = (V_nosecone + V_shoulder) * rho_fiberglass; % [kg]
    end
    
    % -------------- BOATTAIL ---------------
    if strcmp(mode, 'len_bt') || strcmp(mode, 'bt_r2') || i == 1
        len_bt = params.len_bt;
        bt_r2 = params.bt_r2;
        bt_wall_thickness = bt_layers * 0.0003556;
        bt_r1 = ref_diam/2;
        slope_bt = -(bt_r1 - bt_r2) / len_bt;
        
        x2 = linspace(0, len_bt, n_slices);
        y2outer = slope_bt .* x2 + bt_r1;
        y2inner = max(0, y2outer - bt_wall_thickness);
        V_bt = trapz(x2, pi * (y2outer.^2 - y2inner.^2));
        mass_bt = (V_bt + V_shoulder * 4/7) * rho_fiberglass;
        offset_bt = len_bodytubes + len_nose;
        taper_angle = atan((bt_r1-bt_r2) / len_bt);
    end
    % --------------- FINS -----------------
    if strcmp(mode, 'semi_span') || strcmp(mode, 'c') || strcmp(mode, 'tip') || strcmp(mode, 'sweep') || i == 1
        semi_span = params.semi_span;
        c = params.c;
        tip = params.tip;
        sweep_length = params.sweep;
        n_slices2 = 100;
        x3 = linspace(0, semi_span, n_slices2);
        len_chord = ((tip-c) / semi_span) .* x3 + c;
        thick_foam = 0.45 * in_to_m;
        if fintype == 1
            thick_fin = thick_foam+fin_layers*2*thick_cf;
        else
            thick_fin = fin_thickness;
        end
        fin_area = 0.5 * (c + tip) * semi_span;
        fin_area = fin_area + len_bt*(ref_diam/2-bt_r2)/2;
        sweep_angle = atan(sweep_length / semi_span);
        sweep_angle = rad2deg(sweep_angle);
        asp_rat = 2 * semi_span^2 / fin_area;

        
        
        % objective: take the summation of all airfoil cross sections (each section
        % is its own NACA airfoil)
        if fintype == 1
            [volume_1finfoam, ~, ~] = getFinVolumeNACA(c-fin_layers*2*thick_cf, tip-fin_layers*2*thick_cf, thick_foam, semi_span-fin_layers*2*thick_cf);
            [volume_1fin,~,~] = getFinVolumeNACA(c, tip, thick_fin, semi_span);
            mass_foam = volume_1finfoam * rho_foam;
            mass_cf = (volume_1fin-volume_1finfoam) * rho_cf;
            mass_1fin = mass_cf + mass_foam;
            mounting_mass = 4 * (0.05 + 0.025 + 0.005); % mass from fillet, L-bracket and screw respectively
        else
            volume1fin = fin_area*thick_fin*0.85; % just a guess for now FIX THIS LATER FOR VOLUME OF WEDGE AIRFOILS
            mass_1fin = volume1fin * rho_alum;
            mounting_mass = 4 * 0.3;
        end
        surface_area = 2.02 * fin_area; % approximation / slightly over double the planform area
        S_fins_exposed = surface_area * 4; 
        S_fins = S_fins_exposed + (4*c*ref_diam/2); % area for each fin plus the imaginary extension to the center of the body tube
        mass_fins = 4 * mass_1fin;

        % mounting mass
        mass_fins = mass_fins + mounting_mass;
        
        mean_chord = (c + tip) / 2;
        % midslope calc for lift slope
        mid_root = c/2;
        mid_tip = sweep_length + tip/2;
        midslope = (mid_tip - mid_root) / semi_span;
        offset_fins = offset_bt + len_bt - c;
        radius_eff = bt_r2 + tan(taper_angle) * c/2; % effective radius at mid root chord

    end
    rocket_length = len_bodytubes + len_nose + len_bt;
    len_forebody = len_bodytubes + len_nose;
    
    
    
    k = 20e-6; % roughness height in meters
    
    
    x_basedrag = linspace(0, len_bt, n_slices);
    S_bt = trapz(x2, y2outer * 2*pi);
    S_body = len_bodytubes * ref_diam * pi;


    %S_nose = trapz(x, 2 * pi * r_outer); old formula for S_nose
    S_wet = S_nose + S_body + S_bt;
    
    % transonic adjustment estimate: based on Mandell page 480
    M_points = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.4, 1.6];
    transonic_drag = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1.01, 1.05, 1.1, 1.4, 1.7, 1.65, 1.6, 1.57, 1.4, 1.35];
    
    rocket_mass = rocket_mass_no_aero + mass_fins + mass_bt + mass_nosecone;
    
    %% Pack variables into params struct
    params.len_nose = len_nose;
    params.rocket_mass = rocket_mass;
    params.rocket_length = rocket_length;
    params.mean_chord = mean_chord;
    params.k = k;
    params.S_fins = S_fins;
    params.A_ref_nose = A_ref_nose;
    params.thick_fin = thick_fin;
    params.S_fins_exposed = S_fins_exposed;
    params.len_forebody = len_forebody;
    params.dia_nose = ref_diam;
    params.S_wet = S_wet;
    params.bt_r2 = bt_r2;
    params.M_points = M_points;
    params.transonic_drag = transonic_drag;
    params.burn_time = burn_time;
    params.g = g;
    params.S_nose = S_nose;
    params.prop_mass = prop_mass;
    current_val = params.(mode); % essentially a pointer, looks for variable that matches what is inside mode
    params.len_tank = len_tank;
    params.bt_r1 = bt_r1;
    params.c = c;
    params.tip = tip;
    params.sweep = sweep_length;
    params.sweep_angle = sweep_angle;
    params.fin_area = fin_area;
    params.asp_rat = asp_rat;
    params.midslope = midslope;
    params.radius_eff = radius_eff;
    params.mass_fins = mass_fins;
    params.CG_no_fins = CG_no_fins;
    params.offset_bt = offset_bt; 
    params.rasmach = rasmach;
    params.rascd = rascd;
    params.getRasAeroDragImpulse = getRasAeroDragImpulse;
    params.V_wave_haack = V_wave_haack;
    %%
   
    
    % ode45 apogee calculation
    tspan = [0 100];     % Run for 50 seconds
    y0 = [0; 0; 0];        % [initial_alt; initial_vel]
    
    apo = odeset('Events', @apogee_event, 'RelTol', 1e-6, 'AbsTol', 1e-8);    
    [t_out, y_out, te] = ode45(@(t,y) flight(t, y, thrust_time, thrust_force, params), tspan, y0, apo);
    
    % Extract results
    final_alt = y_out(:,1);
    final_vel = y_out(:,2);
    plot(t_out,final_vel) 
    hold on;
    plot(t_out,final_alt)
    rasdragimpulse = y_out(:,3);
    apogees(i) = max(final_alt) * 3.2808399; % apogee in feet
    masses(i) = mass_nosecone / lbs_to_kg; % CHANGE THIS LINE IF YOU WANT TO CHECK NEW MASS
    mass_totals(i) = rocket_mass / lbs_to_kg; 
    if getRasAeroDragImpulse
        rasdragimpulsetotal = max(rasdragimpulse);
        massimpulsenose = mass_bt*g*t_out(end);
        total_impulse = rasdragimpulsetotal + massimpulsenose;
    end
    if plotDamping
        min_zeta_iteration = inf; % Initialize with a high number
    
        for j = 1:length(t_out)
            [rho_j, a_j] = atmo(final_alt(j));
            Q_j = 0.5 * rho_j * final_vel(j)^2;
            if Q_j > max_Q
                max_Q = Q_j; 
                if max_Q > 1.5e5
                    disp(final_vel(j)/a_j)
                end
            end
            % Call stability for each point in the flight
            [current_zeta,~,~] = get_stability(Q_j, rho_j, final_vel(j), t_out(j), params);
            
            % Track the lowest value seen in this specific flight
            if current_zeta < min_zeta_iteration
                min_zeta_iteration = current_zeta;
            end
        end
        min_dampings(i) = min_zeta_iteration;
    end


    if (plotDampingProfile == 1) && (abs(semi_span - target_semi_span) < 1e-5)  
        dampingProfileTime = t_out;
        for j = 1:length(dampingProfileTime)
            [rho_j, a_j] = atmo(final_alt(j));
            Q_j = 0.5 * rho_j * final_vel(j)^2;

            % Call stability for each point in the flight
            [current_zeta,~,~] = get_stability(Q_j, rho_j, final_vel(j), t_out(j), params);
            
            profile_dampings(j) = min(current_zeta,0.15);
        end
        avg_damping = trapz(t_out, profile_dampings)/t_out(end)
    end
    
    % Store the results for plotting

    params.(mode) = current_val + step_size;

end % end of MASSIVE for loop


values = linspace(lower_bound, upper_bound, step_count);
values = values / in_to_m;
% APOGEE PLOT
if plotApogee
    values = linspace(lower_bound, upper_bound, step_count);
    values = values / in_to_m;
    % 
    plot(values, apogees, 'r', 'LineWidth', 2);
    grid on;
    xlabel([mode, '(in)']);
    ylabel('Apogee (ft)');
    title(['Optimizing: ', 'Nosecone length']);
end

% MASS PLOT
if plotCompMass
    figure('WindowStyle','normal');
    plot(values, masses, 'y', 'LineWidth', 2);
    hold on;
    grid on;
    xlabel(mode);
    ylabel('mass(lbm)');
    title(['Mass of : ', mode]);
end

% DAMPING PLOT
if plotDamping
    figure('WindowStyle','normal'); %#ok<*UNRCH>
    plot(values, min_dampings, 'g', 'LineWidth', 2);
    hold on;
    grid on;
    xlabel(mode);
    ylabel('Min Damping Ratio (\zeta)');
    title(['Damping Analysis: ', mode]);
end

% Profile DAMPING PLOT
if plotDampingProfile
    figure('WindowStyle','normal');
    plot(dampingProfileTime, profile_dampings, 'g', 'LineWidth', 2);
    hold on;
    grid on;
    xlabel('time');
    ylabel('Damping Ratio (\zeta)');
    title(['Damping Analysis: ', mode]);
end
%%
function [rho, a] = atmo(alt)
    % Constants
    P0 = 101325;    % Sea level pressure (Pa)
    T0 = 288.15;    % Sea level temp (K) - Standard Atmosphere
    L  = 0.0065;    % Temperature lapse rate (K/m)
    R  = 287.05;    % Gas constant for air (J/kg*K)
    g  = 9.80665;   % Gravity (m/s^2)
    gamma = 1.4;    % Ratio of specific heats

    % Current Temperature (Standard Lapse Model)
    temp = T0 - L * alt;
    temp = max(temp, 216.65); % Floor at Tropopause temp to prevent errors
    a = sqrt(gamma * R * temp);
    pressure = P0 * (1 - L*alt/T0)^(g/(L*R));
    rho = pressure / (R * temp);
end

function dydt = flight(t, h, thrust_time, thrust_force, params)
    %disp(params.S_wet)
    %disp(params.rocket_mass)
    alt = h(1);
    vel = h(2);
    [rho, a] = atmo(alt);
    mach = vel / a;
    
    
    % get interpolated thrust values
    current_thrust = interp1(thrust_time, thrust_force, t, "linear", 0);
    % get current mass

    m_dot = params.prop_mass / params.burn_time;
    if(t < params.burn_time)
        current_mass = (params.rocket_mass - (m_dot*t));
    else
        current_mass = (params.rocket_mass - params.prop_mass);
    end
    % get drag forces
    Re = rho * vel * params.rocket_length / (1.789e-5);
    Re_fin = rho * vel * params.mean_chord / (1.789e-5);
    

    % ~~~~~~~~~~~~~~ Fin drag !!! ~~~~~~~~~~~~~~~~
    % 1. Smooth Turbulent Friction (Depends on Re)
    Cf_smooth = 0.455 / (log10(max(Re, 1e5))^2.58);
    
    % 2. Fully Rough Friction (Depends on Roughness Height k)
    Cf_rough = (1.89 + 1.62 * log10(params.rocket_length / params.k))^-2.5;
    
    % 3. Final Cf
    C_f = max(Cf_smooth, Cf_rough);
    
    % Repeat for Fins
    Cf_fin_smooth = 0.455 / (log10(max(Re_fin, 1e5))^2.58);
    Cf_fin_rough = (1.89 + 1.62 * log10(params.mean_chord / params.k))^-2.5;
    C_f_fins = max(Cf_fin_smooth, Cf_fin_rough);
    C_d_fins = 2 * C_f_fins * (1 + 2 * (0.09)) * (params.S_fins/params.A_ref_nose);
    C_d_interference = 2 * C_f_fins * (1 + 2 * params.thick_fin / params.mean_chord) * ((params.S_fins - params.S_fins_exposed)/params.A_ref_nose);

    % ~~~~~~~~~~~~~~ Body drag !!! ~~~~~~~~~~~~~~~~~~~~~~
    % fineness ratios for nosecone and forebody
    f_n = params.len_nose / params.dia_nose;
    f_b = params.len_forebody / params.dia_nose;
    C_df_body = C_f * (1 + 60/(f_b^3) + 0.0025*f_b) * ((params.S_nose) / params.A_ref_nose);
    C_df_total = C_df_body;

   

    if t > params.burn_time
        C_db = (0.029 / sqrt(C_df_total)) * (params.bt_r2*2 / params.dia_nose)^3;
    else
        C_db = 0;
    end
    if mach > 1.05
        Cd_wave = 9*pi/2 * (params.V_wave_haack/params.len_nose^3)^2 / min(sqrt(max(mach, 1.15)^2-1), 0.5);
        C_df_total = C_df_total + Cd_wave;
    end
    C_d_basic = C_db + C_d_fins + C_d_interference + C_df_total;
    


    



    % ~~~~~~~~~~~~~~~ Corrections for transonic wave drag !!!
    transonic_factor = interp1(params.M_points, params.transonic_drag, mach, 'pchip', 'extrap');
    
    if mach>0.8
        C_d = C_d_basic * transonic_factor;
    else
        C_d = C_d_basic;
    end
    Q_current = 0.5 * rho * vel ^ 2;
    drag_force = Q_current * C_d * params.A_ref_nose;
    
    % Rasaero drag impulse 
    %%
    if params.getRasAeroDragImpulse 
        if mach>0
            targetMach = ceil(mach*100)/100; % discretize current mach to match rasmach
            idx = find(targetMach == params.rasmach, 1, 'first'); % find the index in rasmach that corresponds to current mach
            targetcd = params.rascd(idx);
            rasdragforce = targetcd * Q_current*params.A_ref_nose;
        else 
            rasdragforce = 0;
        end
    else 
        rasdragforce = 0;
    end
    
    % a = F/m
    total_force = current_thrust - (current_mass * params.g) - drag_force;
    acc = total_force / current_mass;
    dydt = [vel; acc; rasdragforce];

    
    
end

function [value, isterminal, direction] = apogee_event(t,h)
    value = h(2);    
    isterminal = 1;    
    direction = -1;
end

function [coupled_damping_ratio, damping_ratio, static_stability] = get_stability(Q, rho, u_vert, time, params)
    global in_to_m;
    
    % Compute real instantaneous Mach number
    temp = max(288.15 - 0.0065 * 0, 216.65); % simplify or pass 'a' from flight
    a = sqrt(1.4 * 287.05 * temp);
    mach = max(0.01, min(0.95, u_vert / a)); % Cap subsonic Prandtl-Glauert below Mach 1
    
    beta = sqrt(1 - mach^2);
    
    FD_term = params.asp_rat / (cos(params.midslope) / beta);
    CN_alpha_alt_num = 2 * 2 * pi * FD_term * (params.fin_area / params.A_ref_nose) * cos(params.midslope);
    CN_alpha_alt_denom = 2 + FD_term * sqrt(1 + (4 / FD_term^2));
    CNA_pre = CN_alpha_alt_num / CN_alpha_alt_denom;
    
    tau = (params.semi_span + params.radius_eff) / params.radius_eff;
    interf_coef = 1 + 1 / tau;
    CN_alpha_fins = CNA_pre * interf_coef;
    
    CN_alpha_nose = 2;
    CN_alpha_boattail = 2 * ((params.bt_r2 / (params.dia_nose/2))^2 - (0.5)^2);
    CN_alpha = CN_alpha_boattail + CN_alpha_nose + CN_alpha_fins;
    
    % Center of Pressure (CP)
    offset_fins = params.offset_bt + params.len_bt - params.c;
    CP_fins = offset_fins + params.sweep/3 * ((params.c + 2*params.tip) / (params.c + params.tip)) + 1/6 *(params.c + params.tip - (params.c*params.tip) / (params.c + params.tip));
    CP_nose = 0.5 * params.len_nose;
    CP_bt = params.offset_bt + (params.len_bt/3) * (1 + params.bt_r2/params.bt_r1 * (params.bt_r2/params.bt_r1 + 1));
    CP_total = (CN_alpha_fins * CP_fins + CN_alpha_nose * CP_nose + CN_alpha_boattail * CP_bt) / CN_alpha;
    
    % Center of Gravity (CG)
    rocket_mass_no_fins = params.rocket_mass - params.mass_fins;
    CG_with_fins = params.CG_no_fins + (params.mass_fins * (offset_fins + 0.54*params.c)) / (params.mass_fins + rocket_mass_no_fins);
    
    m_dot = params.prop_mass / params.burn_time;
    tank_top = params.offset_bt + params.len_bt - params.len_tank; 
    fuel_cg_static = tank_top + params.len_tank / 2;
    
    M_wet = rocket_mass_no_fins + params.mass_fins; 
    m_burned = min(params.prop_mass, m_dot * time);
    cg_dynamic = (M_wet * CG_with_fins - m_burned * fuel_cg_static) / (M_wet - m_burned);
   
    % Aerodynamic derivatives
    C1 = Q * params.A_ref_nose * CN_alpha * (CP_total - cg_dynamic);
    C2A = (rho / 2) * abs(u_vert) * params.A_ref_nose * ( ...
        CN_alpha_fins * (CP_fins - cg_dynamic)^2 + ...
        CN_alpha_nose * (CP_nose - cg_dynamic)^2 + ...
        CN_alpha_boattail * (CP_bt - cg_dynamic)^2 ...
    );
    
    C2B = 0;
    if time <= params.burn_time
        C2B = m_dot * (params.rocket_length - cg_dynamic)^2;
    end
    C2 = C2A + C2B;
    
    % Longitudinal Inertia (IL) including fin pitch inertia
    if time <= params.burn_time
        IL_no_fins = 195 - time * (195-180) / params.burn_time;
    else
        IL_no_fins = 180;
    end
    
    % Correct Fin Pitch Inertia: m * d^2
    d_fin_cg = abs((offset_fins + 0.5 * params.c) - cg_dynamic);
    IL_fins = params.mass_fins * d_fin_cg^2; 
    IL = IL_no_fins + IL_fins;
    
    IR = 0.70; % Roll inertia
    
    % Scalar Damping Ratios
    damping_ratio = C2 / (2 * sqrt(C1 * IL));
    coupled_damping_ratio = C2 / (2 * sqrt(C1 * (IL + IR)));
    static_stability = (CP_total - cg_dynamic) / params.dia_nose;
end



function [finvol, CG_x, avg_thick] = getFinVolumeNACA(c, tip, thick_fin, semi_span)
    n_slices = 100;
    spanwise_areas = zeros(1, n_slices);
    local_centroids = zeros(1, n_slices);
    x1 = linspace(0, semi_span, n_slices);
    len_chord = ((tip - c) / semi_span) .* x1 + c;
    
    offsets = linspace(0, c - tip, n_slices);
    
    for i = 1:n_slices
        x2 = linspace(0, len_chord(i), 100);
        term = x2 ./ len_chord(i);
        airfoil_shape = 5 .* thick_fin .* (0.2969*sqrt(term) - 0.1260*term ...
            - 0.3516*term.^2 + 0.2843*term.^3 - 0.1015*term.^4);
    
        area = 2 * trapz(x2, airfoil_shape);
        spanwise_areas(i) = area;
    
        local_x_cent = trapz(x2, x2 .* airfoil_shape) / trapz(x2, airfoil_shape);
    
        local_centroids(i) = offsets(i) + local_x_cent;
    
    end
    
    
    term = x2 ./ len_chord(i);
    airfoil_shape = 5 .* thick_fin .* (0.2969*sqrt(term) - 0.1260*term ...
        - 0.3516*term.^2 + 0.2843*term.^3 - 0.1015*term.^4);
    avg_thick = 2 * trapz(x2, airfoil_shape) / len_chord(i);
    
    finvol = trapz(x1, spanwise_areas);
    CG_x = trapz(x1, spanwise_areas .* local_centroids) / finvol;
end
% alp = atan(radius_eff/len_nose);
% 1.4 * 2*pi*(nc_layers*thick_fg)^2 * 83e9*cos(alp)^2/sqrt(3*(1-0.21^2))