% wave drag calculator for supersonic airfoils with customization options
% for leading and trailing edge angles and max thickness
% Author: Logan Trickey
% references:
% https://seitzman.gatech.edu/classes/ae2010/supersonicairfoils.pdf
% https://www.me.psu.edu/cimbala/me320/Lesson_Notes/Fluid_Mechanics_Lesson_15F.pdf
% https://www.grc.nasa.gov/www/k-12/airplane/pranmyer.html
% https://www.youtube.com/watch?v=u2oSK19-2ao

% input c (root chord), semi_span, AoA (angle of attack), M_inf
% (free-stream mach number), and p_inf (freestream pressure)
% also enter airfoil dimensions and shape

clear; clc;
gama = 1.4;
in_to_m = 0.0254;

%% ENTER PARAMETERS HERE
c = 18 * in_to_m; 
AoA = deg2rad(10);
semi_span = 7*in_to_m;
M_inf = 2;
p_inf = 80000;
% leading edge and trailing edge **HALF** angles for either airfoil type
angle_LE = deg2rad(5);
angle_TE = deg2rad(5);
fin_thickness = 0.3 * in_to_m;
airfoilshape = 'double wedge';
LE_range = [3.5, 4]; % in degrees
TE_range = [3 5];
DW_min_mid_section = 14*in_to_m; % minimum length for the flat section of double wedge airfoil because as the 
% length of the middle section approaches zero, the math does not accurately treat it as an effective diamond airfoil


step_size = 0.025;
LE_angles = LE_range(1):step_size:LE_range(2);
LE_trial_count = length(LE_angles);


if strcmp(airfoilshape, 'diamond')
    % diamond airfoil geometry
    drag_vec = zeros(1, LE_trial_count);
    TE_angles = zeros(1, LE_trial_count);
    
    for i = 1:LE_trial_count
        angle_LE = deg2rad(LE_angles(i));
        
        % Calculate the exact matching TE angle for a diamond profile
        tan_TE = (fin_thickness * tan(angle_LE)) / (c * tan(angle_LE) - fin_thickness);
        angle_TE = atan(tan_TE);
        TE_angles(i) = rad2deg(angle_TE);
        
        if c - fin_thickness/2 * (1/tan(angle_LE) + 1/tan(angle_TE)) < 0 || tan_TE <= 0
            drag_vec(i) = NaN;
            continue;
        end
        
        drag_vec(i) = airfoildrag(airfoilshape, M_inf, AoA, gama, p_inf, c, fin_thickness, angle_LE, angle_TE, semi_span);
    end
    
    % Plot 2D line for Diamond Airfoil
    figure;
    plot(2*LE_angles, drag_vec, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('Leading Edge Angle (degrees)');
    ylabel('Drag Coefficient');
    title('Diamond Airfoil: Drag Coefficient vs Leading Edge Angle');
    
else
    % double wedge geometry
    TE_angles = TE_range(1):step_size:TE_range(2);
    TE_trial_count = length(TE_angles);
    
    drag_matrix = zeros(LE_trial_count, TE_trial_count);
    
    for i = 1:LE_trial_count
        for j = 1:TE_trial_count
            angle_LE = deg2rad(LE_angles(i));
            angle_TE = deg2rad(TE_angles(j));
            
            if c - fin_thickness/2 * (1/tan(angle_LE) + 1/tan(angle_TE)) < DW_min_mid_section
                drag_matrix(i, j) = NaN;
                continue;
            end
            
            drag_matrix(i, j) = airfoildrag(airfoilshape, M_inf, AoA, gama, p_inf, c, fin_thickness, angle_LE, angle_TE, semi_span);
        end
    end
    
    % Plot 3D Surface for Double Wedge Airfoil
    [TE_grid, LE_grid] = meshgrid(TE_angles, LE_angles);
    
    figure;
    surf(2*LE_grid, 2*TE_grid, drag_matrix);
    xlabel('Leading Edge Angle (degrees)');
    ylabel('Trailing Edge Angle (degrees)');
    zlabel('Drag Coefficient');
    title('Double Wedge: Drag Coefficient Variation with Edge Angles');
    colorbar;
end
%%


function [cd,dragforce] = airfoildrag(airfoilshape, M_inf, AoA, gama, p_inf, c, t, angle_LE, angle_TE, semi_span)

    % Define surface inclination angles relative to freestream flow vector
    if strcmp(airfoilshape, 'double wedge')
        % Rows: 1 = Upper, 2 = Lower
        % Columns: 1 = Front, 2 = Mid, 3 = Rear
        flow_angles = [angle_LE - AoA, -AoA, -angle_TE - AoA; 
                       angle_LE + AoA,  AoA, -angle_TE + AoA];
                   
        % Shoulder turn angles for expansion fans
        turn_angles = [0, angle_LE, angle_TE;
                       0, angle_LE, angle_TE];
                   
        len_mat = [t/(2*sin(angle_LE)), c - t/(2*tan(angle_LE)) - t/(2*tan(angle_TE)), t/(2*sin(angle_TE));
                   t/(2*sin(angle_LE)), c - t/(2*tan(angle_LE)) - t/(2*tan(angle_TE)), t/(2*sin(angle_TE))];
                   
    elseif strcmp(airfoilshape, 'diamond')
        flow_angles = [angle_LE - AoA, -angle_TE - AoA;
                       angle_LE + AoA, -angle_TE + AoA];
                   
        turn_angles = [0, angle_LE + angle_TE;
                       0, angle_LE + angle_TE];
                   
        len_mat = [t/(2*sin(angle_LE)), t/(2*sin(angle_TE));
                   t/(2*sin(angle_LE)), t/(2*sin(angle_TE))];
    end

    P_mat = zeros(size(flow_angles));
    M_mat = zeros(size(flow_angles));

    % State calculations
    for i = 1:size(flow_angles, 1)
        for j = 1:size(flow_angles, 2)
            if j == 1
                % Leading edge: compression shock or expansion fan based on flow angle
                if flow_angles(i,j) > 0
                    [P_mat(i,j), M_mat(i,j)] = oblique(M_inf, flow_angles(i,j), gama, p_inf);
                elseif flow_angles(i,j) < 0
                    [P_mat(i,j), M_mat(i,j)] = expfan(M_inf, abs(flow_angles(i,j)), gama, p_inf);
                else
                    P_mat(i,j) = p_inf;
                    M_mat(i,j) = M_inf;
                end
            else
                % Downstream surfaces use local shoulder turn angles
                [P_mat(i,j), M_mat(i,j)] = expfan(M_mat(i,j-1), turn_angles(i,j), gama, P_mat(i,j-1));
            end
        end
    end

    % Integrate wave drag coefficient (projecting pressure onto freestream vector)
    P_ratio = P_mat ./ p_inf;
    cd = (2 / (gama * M_inf^2 * c)) * sum(sum(P_ratio .* len_mat .* sin(flow_angles)));
    dragforce = cd * 0.5 * gama * p_inf * M_inf^2 * c * semi_span;
end




function [P_abs_oblique, M2] = oblique(M, AoA, gama, P_init)
    theta = abs(AoA);
    [b1,t1] = calc_beta_max(M, gama);
    mu = asin(1/M); % Mach angle
    
    shockAngleFunc = @(B) 2*cot(B) * ((M^2)*sin(B)^2-1) / (M^2*(gama+cos(2*B)) + 2) - tan(theta);
    
    if theta < t1 
        % Bracket weak shock between Mach angle and B_max
        BetaFinal = fzero(shockAngleFunc, [mu + 1e-6, b1 - 1e-6]);
    else 
        error('Aero:ShockDetached', ...
              'The method is invalid due to detached shock. Angle theta (%.2f deg) exceeds theta_max (%.2f deg) at M = %.2f. Try increasing your mach or decreasing your leading edge angle', ...
              rad2deg(theta), rad2deg(t1), M);
    end
    P_rel_oblique = 2*gama/(gama+1) * M^2*sin(BetaFinal)^2-(gama-1)/(gama+1); % pressure gradient from oblique shock

    Mn1 = M*sin(BetaFinal);
    Mn2 = sqrt((1+((gama-1)/2)*Mn1^2) / (gama*Mn1^2-(gama-1)/2));
    M2 = Mn2/sin(BetaFinal-theta);
    P_abs_oblique = P_rel_oblique * P_init;
end

function [P_abs_expansion, M2] = expfan(M, AoA, gama, P_init)
    theta = AoA;
    M2 = prandtlmeyer(M, gama, theta);
    P_rel_expansion = ( (1+(gama-1)/2*M^2) / (1+(gama-1)/2*M2^2) )^ (gama/(gama-1));
    P_abs_expansion = P_rel_expansion * P_init;
end

function [B_max, T_max] = calc_beta_max(M_inf, gama)
    sinesquareterm = 1/(gama*M_inf^2) * ...
       ( (gama+1)/4 * M_inf^2 - 1 + sqrt((gama+1)*(1+(gama-1)/2*M_inf^2+(gama+1)/16*M_inf^4)));
    B_max = asin(sqrt(sinesquareterm));

    T_max = atan(2*cot(B_max) * ((M_inf^2)*sin(B_max)^2-1)  /  (M_inf^2*(gama+cos(2*B_max)) + 2) );
end

function M2 = prandtlmeyer(M1, gama, theta)
    nu = @(M) sqrt((gama+1)/(gama-1))*atan(sqrt((gama-1)/(gama+1)*(M^2-1)))-atan(sqrt(M^2-1));
    nu2_target = nu(M1) + theta;
   
    pm_root_func = @(M2_val) nu(M2_val) - nu2_target;
    M2 = fzero(pm_root_func, [M1, M1 + 10]); % finds the root of the PM function to get the post-expansion mach with initial guess M1 to M1+10
end
