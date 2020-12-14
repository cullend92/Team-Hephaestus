%===================================================================================================
%[]FUNCTION NAME: StateToCoe.m
%[]AUTHOR: Julio César Benavides
%[]CREATED: 06/24/2005
%[]REVISED: 05/30/2016
%===================================================================================================
%[]FUNCTION DESCRIPTION:
%This function calculates the equivalent relative two-body problem classical orbital elements for a
%given set of state vectors.
%===================================================================================================
%[]INPUT VARIABLES:
%(S)|State matrix [km,km/s].
%---------------------------------------------------------------------------------------------------
%(u)|Sum of the gravitational parameters [km^3/s^2].
%===================================================================================================
%[]OUTPUT VARIABLES:
%(coe)|Classical orbital elements [km,-,deg,deg,deg,deg].
%===================================================================================================
%[]VARIABLE FORMAT:
%(S)|Column Vector or Matrix {6 x n}.
%---------------------------------------------------------------------------------------------------
%(u)|Scalar {1 x 1}.
%---------------------------------------------------------------------------------------------------
%(coe)|Column Vector or Matrix {6 x n}.
%===================================================================================================
%[]AUXILIARY FUNCTIONS:
%None.
%===================================================================================================
%[]COMMENTS:
%The classical orbital elements determined by this function are the orbit's semimajor axis,
%eccentricity, inclination, right ascension of the ascending node, argument of periapsis, and true
%anomaly.  All states should be input with respect to an inertial coordinate system.
%===================================================================================================
function Coe = StateToCoe(S,u)
    
    s = size(S);
    %[]Returns the dimensions of the state matrix.
    
    Coe = zeros(s);
    %[]Allocates memory for the classical orbital elements matrix.
    
    I = [1; 0; 0];
    %[]Inertial X-axis vector.
    
    K = [0; 0; 1];
    %[]Inertial Z-axis vector.
    
    for k = 1:s(2);
        
        R = S(1:3,k);
        %[km]Current position vector.
        
        r = norm(R);
        %[km]Magnitude of the current position vector.
        
        V = S(4:6,k);
        %[km/s]Current velocity vector.
        
        v = norm(V);
        %[km/s]Magnitude of the current velocity vector.
        
        H = cross(R,V);
        %[km^2/s]Current specific angular momentum vector.
        
        h = norm(H);
        %[km^2/s]Magnitude of the current specific angular momentum.
        
        N = cross(K,H);
        %[km^2/s]Current ascending node vector.
        
        n = norm(N);
        %[m^2/s]Current ascending node magnitude.
        
        Dot = dot(R,V);
        %[km^2/s]Dot product of the current position and velocity vectors.
        
        E = 1 / u * ((v^2 - u / r) * R - Dot * V);
        %[]Current eccentricity vector.
        
        e = norm(E);
        %[]Current eccentricity.
        
        epsilon = v^2 / 2 - u / r;
        %[km^2/s^2]Current specific mechanical energy.
        
        Coe(1,k) = -u / (2 * epsilon);
        %[km]Current semimajor axis.
        
        Coe(2,k) = e;
        %[]Current eccentricity magnitude.
        
        Coe(3,k) = acosd(dot(K,H) / h);
        %[deg]Current inclination.
        
        Coe(4,k) = acosd(dot(I,N) / n);
        %[deg]Current right ascension of the ascending node.
        
        if N(2) < 0;
            
            Coe(4,k) = 360 - Coe(4,k);
            %[deg]Right ascension of the ascending node quadrant check.
            
        end
        
        Coe(5,k) = acosd(dot(N,E) / (n * e));
        %[deg]Current argument of periapsis.
        
        if E(3) < 0;
            
            Coe(5,k) = 360 - Coe(5,k);
            %[deg]Argument of periapsis quadrant check.
            
        end
        
        Coe(6,k) = acosd(dot(E,R) / (e * r));
        %[deg]Current true anomaly.
        
        if Dot < 0;
            
            Coe(6,k) = 360 - Coe(6,k);
            %[deg]True anomaly angle quadrant check.
            
        end
        
    end
    
end
%===================================================================================================