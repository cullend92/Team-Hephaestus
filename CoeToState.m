%===================================================================================================
%[]FUNCTION NAME: CoeToState.m
%[]AUTHOR: Julio César Benavides
%[]CREATED: 06/24/2005
%[]REVISED: 05/30/2016
%===================================================================================================
%[]FUNCTION DESCRIPTION:
%This function calculates the state vectors for an orbit's set of relative two-body problem
%classical orbital elements.
%===================================================================================================
%[]INPUT VARIABLES:
%(coe)|Classical orbital elements matrix [km,-,deg,deg,deg,deg].
%---------------------------------------------------------------------------------------------------
%(u)|Sum of the gravitational parameters [km^3/s^2].
%===================================================================================================
%[]OUTPUT VARIABLES:
%(S)|State vectors [km,km/s].
%===================================================================================================
%[]VARIABLE FORMAT:
%(coe)|Column Vector or Matrix {6 x n}.
%---------------------------------------------------------------------------------------------------
%(u)|Scalar.
%---------------------------------------------------------------------------------------------------
%(S)|Column Vector or Matrix {6 x n}.
%===================================================================================================
%[]AUXILIARY FUNCTIONS:
%None.
%===================================================================================================
%[]COMMENTS:
%Every column vector in the orbital elements matrix should have the six (6) classical orbital
%elements of the body being analyzed (semimajor axis, eccentricity, inclination, right ascension of
%ascending node, argument of perigee, and true anomaly).  Inclination, right ascension of the
%ascending node, argument of perigee, and true anomaly should be in degrees.  'PQW' stands for
%"Perifocal Inertial Coordinate System". 'IJK' stands for "Inertial Coordinate System".  'COE'
%stands for "classical orbital elements". 'WRT' stands for "with respect to".
%===================================================================================================
function S = CoeToState(coe,u)
    
    s = size(coe);
    %[]Returns the dimensions of the COE matrix.
    
    S = zeros(s);
    %[]Allocates memory for the state vectors corresponding the COE matrix.
    
    for k = 1:s(2)
        
        a = coe(1,k);
        %[km]Current semimajor axis.
        
        e = coe(2,k);
        %[]Current eccentricity.
        
        In = coe(3,k) * pi / 180;
        %[rad]Current inclination.
        
        Om = coe(4,k) * pi / 180;
        %[rad]Current right ascension of the ascending node.
        
        w = coe(5,k) * pi / 180;
        %[rad]Current argument of periapsis.
        
        theta = coe(6,k) * pi / 180;
        %[rad]Current true anomaly.
        
        p = a * (1 - e^2);
        %[km]Current semiparameter WRT the central body.
        
        r = p / (1 + e * cos(theta));
        %[km]Current range WRT the central body.
        
        Rpqw = r * [cos(theta); sin(theta); 0];
        %[km]Current position vector WRT the central body in the PQW.
        
        Vpqw = sqrt(u / p) * [-sin(theta); cos(theta) + e; 0];
        %[km/s]Current velocity vector WRT the central body in the PQW.
        
        Spqw = [Rpqw; Vpqw];
        %[km,km/s]Current state vector WRT the central body in the PQW.
        
        R3Om = [cos(Om), -sin(Om), 0; ...
                sin(Om),  cos(Om), 0; ...
                      0,        0, 1];
        %[]Rotation matrix about the third axis by an angle right ascension of the ascending node.
        
        R1In = [1,       0,        0; ...
                0, cos(In), -sin(In); ...
                0, sin(In),  cos(In)];
        %[]Rotation matrix about the first axis by an angle inclination
        
        R3w = [cos(w), -sin(w), 0; ...
               sin(w),  cos(w), 0; ...
                    0,       0, 1];
        %[]Rotation matrix about the third axis by an angle argument of periapsis.
        
        T = R3Om * R1In * R3w;
        %[]Current matrix that transforms vectors from the PQW to the IJK.
        
        S(:,k) = [T, zeros(3,3); zeros(3,3), T] * Spqw;
        %[km,km/s]Current state vector WRT the central body in the IJK.
        
    end
    
end
%===================================================================================================