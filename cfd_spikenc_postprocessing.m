% cfd_drag_calc.m
% Computes drag force over a flight using CFD-derived Cd(Mach) data and
% an OpenRocket (.ork) simulation export containing time, velocity,
% air density, and Mach number.
%
% Workflow:
%   1. Load CFD results: one Cd value per Mach number you simulated.
%   2. Load the .ork export (time, velocity, air density, Mach).
%   3. Interpolate Cd at every timestep's Mach number using the CFD curve.
%   4. Drag(t) = Cd_interp(t) * 0.5 * rho(t) * v(t)^2 * A_ref

clear; clc;

%% ~~~~~~ USER INPUT ~~~~~~

% CFD results: two columns, Mach and Cd, one row per simulated point.
% Edit this path once you have CFD data, or replace with a readmatrix()
% call pointing at your results file.
cfd_file = "cfd_results.csv";      % columns: Mach, Cd (no header)

% OpenRocket export: columns are time [s], total velocity [m/s],
% air density [kg/m^3], Mach (matches the example .ork export, no header)
ork_file = "example ork file containing time totalvelocity airdensity and mach.csv";

% Reference area used to convert Cd into a drag force [m^2].
% Set this to whatever area your CFD Cd values were normalized against
% (usually the rocket's max cross-sectional area).
in_to_m = 0.0254;
rocket_diam = 8.4 * in_to_m;       % rocket diameter [m]
A_ref = pi * (rocket_diam/2)^2;    % max cross-sectional area [m^2]

%% ~~~~~~ LOAD DATA ~~~~~~

cfd_data = readmatrix(cfd_file);
cfd_mach = cfd_data(:,1);
cfd_cd   = cfd_data(:,2);

% interp1 requires strictly increasing x — sort in case CFD runs weren't
% simulated in Mach order
[cfd_mach, sort_idx] = sort(cfd_mach);
cfd_cd = cfd_cd(sort_idx);

% drop duplicate Mach points (interp1 errors on repeated x values) —
% keeps the last Cd entered for any repeated Mach
[cfd_mach, unique_idx] = unique(cfd_mach, 'last');
cfd_cd = cfd_cd(unique_idx);

if numel(cfd_mach) < 2
    error('Need at least 2 distinct CFD Mach points to interpolate — got %d.', numel(cfd_mach));
end

ork_data = readmatrix(ork_file);
time     = ork_data(:,1);
velocity = ork_data(:,2);
rho      = ork_data(:,3);
mach     = ork_data(:,4);

%% ~~~~~~ INTERPOLATE Cd AT EACH TIMESTEP ~~~~~~

mach_min = min(cfd_mach);
mach_max = max(cfd_mach);

% only interpolate — never extrapolate. Timesteps whose flight Mach falls
% outside the simulated CFD range [mach_min, mach_max] get NaN and are
% excluded from the drag calculation entirely.
in_range = mach >= mach_min & mach <= mach_max;

cd_interp = nan(size(mach));
cd_interp(in_range) = interp1(cfd_mach, cfd_cd, mach(in_range), 'pchip');

if any(~in_range)
    fprintf('%d of %d timesteps have Mach outside the CFD range [%.3f, %.3f] — excluded from drag calc\n', ...
        sum(~in_range), numel(mach), mach_min, mach_max);
end

%% ~~~~~~ DRAG FORCE (only within the CFD-simulated Mach range) ~~~~~~

q = 0.5 * rho .* velocity.^2;      % dynamic pressure [Pa]
drag_force = nan(size(mach));
drag_force(in_range) = cd_interp(in_range) .* q(in_range) * A_ref;   % [N]

% trimmed vectors covering only the valid Mach window, for plotting/output
time_v     = time(in_range);
mach_v     = mach(in_range);
cd_v       = cd_interp(in_range);
rho_v      = rho(in_range);
velocity_v = velocity(in_range);
drag_v     = drag_force(in_range);

%% ~~~~~~ PLOTS ~~~~~~

figure;
plot(cfd_mach, cfd_cd, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(mach_v, cd_v, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Mach');
ylabel('C_d');
legend('CFD points', 'Interpolated (flight)', 'Location', 'best');
title('C_d vs Mach');

figure;
plot(time_v, drag_v, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Drag Force (N)');
title(sprintf('Drag Force vs Time (Mach %.2f - %.2f only)', mach_min, mach_max));

%% ~~~~~~ DRAG IMPULSE (within the CFD-simulated Mach range) ~~~~~~

% The flight can pass through the CFD Mach window more than once (e.g.
% accelerating up through it, then falling back through it after
% apogee), so `in_range` may consist of several disjoint time segments.
% Integrating drag_v against time_v directly with trapz would wrongly
% bridge the gap between segments, so each contiguous run is integrated
% separately and the impulses are summed.
idx_in_range = find(in_range);
gap_after = [find(diff(idx_in_range) > 1); numel(idx_in_range)];
seg_start = 1;
segment_impulses = zeros(numel(gap_after), 1);
segment_t_start  = zeros(numel(gap_after), 1);
segment_t_end    = zeros(numel(gap_after), 1);

for k = 1:numel(gap_after)
    seg_idx = idx_in_range(seg_start:gap_after(k));
    segment_t_start(k) = time(seg_idx(1));
    segment_t_end(k)   = time(seg_idx(end));
    if numel(seg_idx) >= 2
        segment_impulses(k) = trapz(time(seg_idx), drag_force(seg_idx));
    else
        segment_impulses(k) = 0;   % can't integrate a single point
    end
    seg_start = gap_after(k) + 1;
end

drag_impulse_total = sum(segment_impulses);   % [N*s]

fprintf('\nDrag impulse within Mach [%.3f, %.3f]: %.4f N*s\n', mach_min, mach_max, drag_impulse_total);
if numel(segment_impulses) > 1
    fprintf('  (split across %d segments where flight re-entered the CFD Mach range)\n', numel(segment_impulses));
    for k = 1:numel(segment_impulses)
        fprintf('    segment %d: t = %.3f - %.3f s, impulse = %.4f N*s\n', ...
            k, segment_t_start(k), segment_t_end(k), segment_impulses(k));
    end
end

%% ~~~~~~ OUTPUT TABLE (Mach range covered by CFD only) ~~~~~~

results = table(time_v, mach_v, cd_v, rho_v, velocity_v, drag_v, ...
    'VariableNames', {'time_s', 'mach', 'Cd', 'rho_kgm3', 'velocity_ms', 'drag_N'});
% writetable(results, 'drag_results.csv');   % uncomment to save
