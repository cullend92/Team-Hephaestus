%   EarthMarsPorkchp.m
%   Authors:
%       Daniel Cullen
%       Aaron Scott
%
%   EP 394 - Spacecraft System Engineering
%   Team Hephaestus - Robotic Assets
%
%   Description:
%       Generates porkchop plots for C3 and arrival delta-V for Earth-mars
%       transfers.
%
%% Preprocessing:
% Parameter Entry:
clear; clc; close all;

c = getPlanetParameters();

% Initial parking orbit parameters: (Circular orbit in LEO)
he = 350; % [km] Initial parking orbit altitude
rpe = he + c.R_earth; % [km] Initial parking orbit perigee
vpi_e = sqrt(c.mu_earth/rpe);

% Target orbit parameters:
e = 0.8; % [] Target orbit eccentricity
h = 315; % [km] Taret orbit altitude
rpm = c.R_mars + h; % [km] Periapsis height of target orbit
a = (rpm)/(1 - e); % [km] Semi-major axis of target orbit
P = a*(1-e^2); % [km] Semiparameter of target orbit
vpf_m = (1 + e)*sqrt(c.mu_mars/P); % [km/s] Periapsis velocity of target orbit

% Trajectory Generation parameters:
numDataPts = 2500;
departureRange_start = [2030, 01, 01]; % [Yr, Mo, D] Departure range start
departureRange_end = [2034, 01, 01]; % [Yr, Mo, D] Departure range end
TOF_min = 90; % [days] Minimum allowable time of flight
TOF_max = 450; % [days] Maximum allowable time of flight
arrivalRange_start = departureRange_start + [0, 0, TOF_min]; % [Yr, Mo, D] Departure range start
arrivalRange_end = departureRange_end + [0, 0, TOF_max];

% Plotting parameters:
dV_max = 15; % [km/s] Max delta V to display
C3_max = 15; % [km^2/s^2] Max C3 to display

% Prepare starting data:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting departureRange state initialization...');
departureRange = linspace(juliandate(departureRange_start),...
                          juliandate(departureRange_end),...
                          numDataPts);
[r1, v1] = planetEphemeris(departureRange', 'Sun', 'Earth', '432t', 'km');
S1 = [r1, v1];
disp('Finished departureRange state initialization.');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Starting arrivalRange state initialization...');
arrivalRange = linspace(juliandate(arrivalRange_start),...
                        juliandate(arrivalRange_end),...
                        numDataPts);
[r2, v2] = planetEphemeris(arrivalRange', 'Sun', 'Mars', '432t', 'km');
S2 = [r2, v2];
disp('Finished arrivalRange state initialization.');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Trajectories:
disp('Calculating trajectories...');
C3 = zeros(numDataPts);
dV_scalar = zeros(numDataPts);
for x = 1:length(departureRange)
    for y = 1:length(arrivalRange)
           TOF = arrivalRange(y) - departureRange(x);
           if (TOF <= 0)
               dV_scalar(y, x) = NaN;
               C3(y, x) = NaN;
           else
               [v_sc_s_depart, v_sc_s_arrive] = ...
                    glambert(c.mu_sun, [r1(x, :) v1(x, :)], [r2(y, :) v2(y, :)], TOF*86400, 0);
                
               vinf_e = v_sc_s_depart - S1(x, 4:6);
               vpf_e = sqrt(norm(vinf_e)^2 + 2*c.mu_earth/rpe);
               dV_depart = vpf_e - vpi_e;
               
               vinf_m = S2(y, 4:6) - v_sc_s_arrive;
               vpi_m = sqrt(norm(vinf_m)^2 + 2*c.mu_mars/rpm);
               dV_arrive = (vpi_m - vpf_m);
               
               C3(y, x) = norm(vinf_e)^2;
               dV_scalar(y, x) = norm(dV_arrive) + norm(dV_depart);
           end
    end
    if mod(x, 10) == 0
        dispProgress(x, numDataPts);
    end
end
disp('Finished trajectory calculation.');

%% Postprocessing:
disp('Beginning Postprocessing.');
% Trim C3 and dV matrices:
C3_display = C3;
C3_mask = (C3 < 0) | (C3 > C3_max);
C3_display(C3_mask) = NaN;

dV_scalar_display = dV_scalar;
dV_mask = (dV_scalar < 0) | (dV_scalar > dV_max);
dV_scalar_display(dV_mask) = NaN;

%% Create delta V Contour plot (porkchop plot):
disp('Generating Figures:');
disp('Generating delta V porkchop plot.');
f1 = figure(1);
contourf(departureRange,...
         arrivalRange,...
         dV_scalar_display,...
         'LineStyle', 'none',...
         'HandleVisibility', 'off');
grid on; hold on;
CB = colorbar('EastOutside');
CB.Label.String = 'km/s';
title('Earth-Mars Transfers Total Delta-V');
xlabel('Departure Date');
ylabel('Arrival Date');

plot(juliandate([2033, 04, 01]), juliandate([2033, 09, 28]), 'g*', 'markersize', 5);
plot(juliandate([2033, 04, 10]), juliandate([2033, 10, 07]), 'r*', 'markersize', 5);
legend('Beginning of departure window', 'End of departure window',...
       'Location', 'SouthEast');

xticks = linspace(departureRange(1), departureRange(end), 15);
xTickLabels = datestr(datetime(xticks, 'convertfrom', 'juliandate'), 1);
yticks = linspace(arrivalRange(1), arrivalRange(end), 15);
yTickLabels = datestr(datetime(yticks, 'convertfrom', 'juliandate'), 1);

set(gca, 'XTick', xticks, 'XTickLabel', xTickLabels);
set(gca, 'YTick', yticks, 'YTickLabel', yTickLabels);
xtickangle(45);

%% Create C3 Contour plot:
disp('Generating C3 porkchop plot.');
f2 = figure(2);
contourf(departureRange,...
         arrivalRange,...
         C3_display,...
         'LineStyle', 'none',...
         'HandleVisibility', 'off');
CB = colorbar('EastOutside');
CB.Label.String = 'km^2/s^2';
title('Earth-Mars Transfers C3');
xlabel('Departure Date');
ylabel('Arrival Date');

set(gca, 'XTick', xticks, 'XTickLabel', xTickLabels);
set(gca, 'YTick', yticks, 'YTickLabel', yTickLabels);
xtickangle(45);

grid on; hold on;
plot(juliandate([2033, 04, 01]), juliandate([2033, 09, 28]), 'g*', 'markersize', 5);
plot(juliandate([2033, 04, 10]), juliandate([2033, 10, 07]), 'r*', 'markersize', 5);
legend('Beginning of departure window', 'End of departure window',...
       'Location', 'SouthEast');
   
disp('Script finished.');

%% Helper function: Display trajectory generation progress
function dispProgress(i_curr, i_max)
    persistent bckspLen
    if isempty(bckspLen)
        bckspLen = 0;
    end
    fprintf(repmat('\b', 1, bckspLen));
    percentComplete = i_curr/i_max*100;
    printString = sprintf(['Trajectory calculations are ',...
                           '%.2f %%%% complete...'],...
                           percentComplete);
    fprintf([printString, '\n']);
    bckspLen = length(printString);
end