clear; clc; close all;

pr = 255; % primary color value
sc = 100; % secondary color values
red = [pr, sc, sc] / norm([pr, sc, sc]);
green = [sc, pr, sc] / norm([sc, pr, sc]);
blue = [sc, sc, pr] / norm([sc, sc, pr]);

% Data Entry:
numDataPts = 1000;
TOF = 180; % [days]
simStart = juliandate([2033, 04, 01]);
simEnd = simStart + TOF;
simRange = linspace(simStart, simEnd, numDataPts);
c = getPlanetParameters();

% Get planet states for use in glambert
[r_earth_launch, v_earth_launch] = planetEphemeris(simStart, 'Sun', 'Earth', '432t', 'km');
[r_mars_arrival, v_mars_arrival] = planetEphemeris(simEnd, 'Sun', 'Mars', '432t', 'km');

% Get planet state matrices for plotting orbit lines
orbitRange = linspace(simStart, simStart + 365, numDataPts);
r_earth = planetEphemeris(orbitRange', 'Sun', 'Earth', '432t', 'km');
orbitRange = linspace(simStart, simStart + 687, numDataPts);
r_mars = planetEphemeris(orbitRange', 'Sun', 'Mars', '432t', 'km');

TOF = TOF * 86400;
[v1, v2] = glambert(c.mu_sun, [r_earth_launch, v_earth_launch],...
                              [r_mars_arrival, v_mars_arrival], TOF, 0);

% Get spacecraft state for plotting trajectory line and animating
% trajectory
coe_sc_start = StateToCoe([r_earth_launch'; v1'], c.mu_sun);
coe_sc_end = StateToCoe([r_mars_arrival'; v2'], c.mu_sun);
coe_array = ones(6, numDataPts) .* coe_sc_start;
coe_array(6, :) = linspace(coe_sc_start(6), coe_sc_end(6), numDataPts);
S_sc = CoeToState(coe_array, c.mu_sun);
S_sc = S_sc';

% Plot sun and earth/mars orbits and trajectory line:
f = figure(1);
set(f, 'units', 'normalized', 'position', [0.01 0.01 0.8 0.8]);
plot3(0, 0, 0, 'y*', 'markersize', 30,...
      'HandleVisibility', 'Off');
hold on; grid on;
plot3(c.AU_per_km * S_sc(:, 1),...
      c.AU_per_km * S_sc(:, 2),...
      c.AU_per_km * S_sc(:, 3), 'color', green,...
      'HandleVisibility', 'Off');
plot3(c.AU_per_km * r_earth(:, 1),...
      c.AU_per_km * r_earth(:, 2),...
      c.AU_per_km * r_earth(:, 3), 'color', blue,...
      'HandleVisibility', 'Off');
plot3(c.AU_per_km * r_mars(:, 1),...
      c.AU_per_km * r_mars(:, 2),...
      c.AU_per_km * r_mars(:, 3), 'color', red,...
      'HandleVisibility', 'Off');
xlim([-2, 2]);
ylim([-2, 2]);
zlim([-2, 2]);
xlabel('X (AU)');
ylabel('Y (AU)');
zlabel('Z (AU)');
set(gca, 'FontSize', 18, 'GridColor', 'w');
title('Earth-Mars Transfer', 'FontSize', 24);

% Generate planet positions for animation:
r_earth = planetEphemeris(simRange', 'Sun', 'Earth', '432t', 'km');
r_mars = planetEphemeris(simRange', 'Sun', 'Mars', '432t', 'km');

r_earth = r_earth * c.AU_per_km;
r_mars = r_mars * c.AU_per_km;
S_sc = S_sc * c.AU_per_km;

EP = plot3(r_earth(1, 1),...
           r_earth(1, 2),...
           r_earth(1, 3), 'bo', 'markersize', 24);
MP = plot3(r_mars(1, 1),...
           r_mars(1, 2),...
           r_mars(1, 3), 'ro', 'markersize', 18);
TVP = plot3(S_sc(1, 1),...
            S_sc(1, 2),...
            S_sc(1, 3), 'go', 'markersize', 12);
dateText = sprintf('Current Date: %s',...
    datestr(datetime(simRange(1), 'convertfrom', 'juliandate'), 1));
T = text(0.75, 0, 0, dateText, 'units', 'normalized', 'FontSize', 24);
legend('Earth', 'Mars', 'Transit Vehicle', 'location', 'northeast', 'FontSize', 18);
set(gca, 'color', 'k');
EP.Color = blue;
MP.Color = red;
TVP.Color = green;

%% Animation
while(1)
    pause(5);
    for i = 1:numDataPts
        TVP.XData = S_sc(i, 1);
        TVP.YData = S_sc(i, 2);
        TVP.ZData = S_sc(i, 3);
        EP.XData = r_earth(i, 1);
        EP.YData = r_earth(i, 2);
        EP.ZData = r_earth(i, 3);
        MP.XData = r_mars(i, 1);
        MP.YData = r_mars(i, 2);
        MP.ZData = r_mars(i, 3);
        dateText = sprintf('Current Date: %s',...
        datestr(datetime(simRange(i), 'convertfrom', 'juliandate'), 1));
        T.String = dateText;
        T.Position = [0.725, 0, 0];
        pause(0.005);
    end
end