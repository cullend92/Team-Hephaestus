%   Retrieves planetary parameters based on values provided by Dr. Kaela
%   Martin.
function [ c ] = getPlanetParameters()
% Gravitational paramater, km^3/sec^2
c.mu_sun = 132712440000;
c.mu_moon = 4902.8;
c.mu_mercury = 22032;
c.mu_venus = 324860;
c.mu_earth = 398600.44;
c.mu_mars = 42828;
c.mu_jupiter = 126713000;
c.mu_saturn = 37941000;
c.mu_uranus = 5794500;
c.mu_neptune = 6836500;
c.mu_pluto = 981.60;
c.mu_titan = 8978;

% Radius, km
c.R_sun = 695990;
c.R_moon = 1738.2;
c.R_mercury = 2439.7;
c.R_venus = 6051.9;
c.R_earth = 6378.136;
c.R_mars = 3397;
c.R_jupiter = 71492;
c.R_saturn = 60268;
c.R_uranus = 25559;
c.R_neptune = 25269;
c.R_pluto = 1162;
c.R_titan = 2574.7;

% Semi-major axis, km
c.a_moon = 384400;
c.a_mercury = 57910000;
c.a_venus = 108210000;
c.a_earth = 149600000;
c.a_mars = 227920000;
c.a_jupiter = 778570000;
c.a_saturn = 1433530000;
c.a_uranus = 2872460000;
c.a_neptune = 4495060000;
c.a_pluto = 5906380000;
c.a_titan = 20.3*c.R_saturn;

% Period, sec
c.P_moon = 235872;
c.P_mercury = 7603200;
c.P_venus = 19414080;
c.P_earth = 31553280;
c.P_mars = 59356800;
c.P_jupiter = 382838400;
c.P_saturn = 928540800;
c.P_uranus = 2642889600;
c.P_neptune = 5166720000;
c.P_pluto = 7826803200;

% Eccentricity, unitless
c.e_moon = 0.0549;
c.e_mercury = 0.2056;
c.e_venus = 0.0067;
c.e_earth = 0.0167;
c.e_mars = 0.0935;
c.e_jupiter = 0.0489;
c.e_saturn = 0.0565;
c.e_uranus = 0.0457;
c.e_neptune = 0.0113;
c.e_pluto = 0.2488;
c.e_titan = 0.0287;

% Inclination w.r.t ecliptic, degrees
c.i_moon = 5.145;
c.i_mercury = 7;
c.i_venus = 3.39;
c.i_earth = 0;
c.i_mars = 1.850;
c.i_jupiter = 1.304;
c.i_saturn = 2.485;
c.i_uranus = 0.772;
c.i_neptune = 1.769;
c.i_pluto = 17.16;

% Unit conversion factors
c.AU_per_km = (1/1.496e8);
c.km_per_AU = (1/6.68459e-9);
end

