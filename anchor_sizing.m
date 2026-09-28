%%% Enter here
tip_length = 1.1438; % in, tip/aerospike length (moment arm from shoulder)
P = 695;             % lbf, max load on the tip/aerospike
d = [0.5 0.75 1 1.25 1.5]; % in, anchor diameters to check
%%%

% material / sheet constants
Sy = 73000;          % psi, 7075 yield (governs over ultimate at the same FOS)
FOS = 1.5;
S_thick = 10000;     % psi, thickened compressive strength

M = P * tip_length;  % lbf*in, moment at shoulder

% shoulder stresses and MOS (yield)
A = pi * (d/2).^2;             % in^2
I = (pi/4) * (d/2).^4;         % in^4
tau = P ./ A;                  % psi, shear stress
sigma = M * (d/2) ./ I;        % psi, bending stress
MOS_shear = Sy ./ (FOS * tau) - 1;
MOS_bend = Sy ./ (FOS * sigma) - 1;
pass = MOS_shear >= 0 & MOS_bend >= 0;

% minimum embedded length (bearing on thickened)
% solves P/(dL)(4+6a/L) = sigma_a
sigma_a = S_thick / FOS;
L_min = (4*P + sqrt(16*P^2 + 24*sigma_a*d*P*tip_length)) ./ (2*sigma_a*d);

% anchor volume at minimum length
V = A .* L_min;                % in^3

fprintf('Moment at shoulder: %.1f lbf*in\n\n', M);
fprintf('  d (in)   Shear MOS   Bend MOS   MOS     L_min (in)   Volume (in^3)\n');
status = {'FAIL', 'pass'};
for i = 1:length(d)
    fprintf('  %5.2f    %8.2f    %7.2f    %-4s    %7.3f       %7.3f\n', ...
        d(i), MOS_shear(i), MOS_bend(i), status{pass(i) + 1}, L_min(i), V(i));
end
