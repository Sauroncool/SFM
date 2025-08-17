% input to the function are the components of initial position and velocity vector
% output are the orbit parameters, that is, eccentricity (e), angular
% momentum of the orbit (h) and true anomaly (thetha0)

function [e,h,thetha0] = orbital_parameters(r0x,r0y,r0z,v0x,v0y,v0z) 

r0 = sqrt((r0x*r0x) + (r0y*r0y) + (r0z*r0z)); % initial position vector magnitude
v0 = sqrt((v0x*v0x) + (v0y*v0y) + (v0z*v0z)); % initial velocity vector magnitude
vr0 = ((r0x*v0x) + (r0y*v0y) + (r0z*v0z))/r0; % radial component of velocity vector magnitude
vp0 = sqrt((v0*v0) - (vr0*vr0)); % perpendicular component of velocity vector magnitude
v0 = sqrt((v0x*v0x) + (v0y*v0y) + (v0z*v0z));
mu = 398600; 

A = (((v0*v0) - (vr0*vr0))*r0/mu) - 1; % e*cos(thetha0) value
B = (vr0*r0/mu)*sqrt(((v0*v0) - (vr0*vr0))); % e*sin(thetha0) value

e = sqrt((A*A) + (B*B));
tt = acosd(A/e); % true anomaly

% Above we will get two values of tt. Now which one to accept depends upon
% the flight path angle, which is calculated as shown below.
gamma = atand(vr0/vp0);
if gamma > 0
    if A > 0
        thetha0 = tt;
        fprintf("1st quad, t = %d", thetha0);
    elseif A < 0
        thetha0 = tt;
        fprintf("2nd quad, t = %d", thetha0);
    end
elseif gamma < 0
    if A > 0
        thetha0 = (-360+tt);
        fprintf("4th quad, t = %d", thetha0);
    elseif A < 0
        thetha0 = (-360+tt);
        fprintf("3rd quad, t = %d", thetha0);
    end
end

h = (mu*e*sind(thetha0))/vr0; % angular momentum calculation

end