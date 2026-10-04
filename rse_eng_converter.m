clear; clc;
str = fileread("C:\Users\ltric\Downloads\Thrust_Curve (1).RSE");

% Find all matches in the giant string
t_matches = regexp(str, ' t="([\d.]+)"', 'tokens');
% Added [eE\-\+]* to allow for scientific notation like 1.23e-05
f_matches = regexp(str, 'f="([\d.eE\-\+]+)"', 'tokens');
%name_match = regexp(str, '')

% Convert from cell array of cells to a simple double array
time = cellfun(@(x) str2double(x{1}), t_matches);
thrust = cellfun(@(x) str2double(x{1}), f_matches);

% Plot to verify the curve
plot(time, thrust);
xlabel('Time (s)'); ylabel('Thrust (N)');
title('Hybrid Motor Thrust Curve');

len = length(time); 

engfile = strings(len, 1); 

for i = 1:len
    engfile(i) = fprintf('%.4f %.4f\n', time(i), thrust(i));
end

% To see the result:
disp(engfile);