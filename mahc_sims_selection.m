clear; clc;

% load flight data (time, velocity, air density, Mach)
data = readmatrix(fullfile(fileparts(mfilename('fullpath')), '..', 'excel misc', 'flight ork time vs vel vs rho vs mach.csv'));
time = data(:, 1);
vel  = data(:, 2);
rho  = data(:, 3);
mach = data(:, 4);

% cut off at apogee (first point after max Mach where Mach starts rising again)
[~, iMax] = max(mach);
iApogee = iMax + find(diff(mach(iMax:end)) > 0, 1) - 1;
time = time(1:iApogee);
vel  = vel(1:iApogee);
rho  = rho(1:iApogee);
mach = mach(1:iApogee);

% weights: time spent in each timestep, and drag (rho * v^2 * dt) in each timestep
dt   = [diff(time); 0];
drag = rho .* vel.^2 .* dt;

% method 1: 8 points spread by drag
disp('Drag weighted:')
disp(pick(mach, drag, 8))

% method 2: 8 points spread by time
disp('Time weighted:')
disp(pick(mach, dt, 8))

% split flight into subsonic, transonic, supersonic
sub   = mach < 0.8;
trans = mach >= 0.8 & mach <= 1.2;
sup   = mach > 1.2;

% method 3: 3 subsonic, 3 transonic, 2 supersonic, spread by drag within each
disp('Regime drag weighted:')
disp([pick(mach(sub),   drag(sub),   3), ...
      pick(mach(trans), drag(trans), 3), ...
      pick(mach(sup),   drag(sup),   2)])

% method 4: 3 subsonic, 3 transonic, 2 supersonic, spread by time within each
disp('Regime time weighted:')
disp([pick(mach(sub),   dt(sub),   3), ...
      pick(mach(trans), dt(trans), 3), ...
      pick(mach(sup),   dt(sup),   2)])

% picks n Mach numbers, each in the middle of an equal share of the weight
function machCFD = pick(mach, weight, n)
    % sort by Mach and build running total of weight (0 to 1)
    [machSorted, order] = sort(mach);
    total = cumsum(weight(order));
    total = total / total(end);

    % middle of each of the n equal slices, then the Mach where the total reaches it
    targets = ((1:n) - 0.5) / n;
    machCFD = zeros(1, n);
    for k = 1:n
        machCFD(k) = machSorted(find(total >= targets(k), 1));
    end
end
