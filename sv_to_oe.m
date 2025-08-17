function [h, e, ra, inc, w, ta, a] = sv_to_oe(R,V)

mu = 398600;
r = norm(R); % magnitute of radius vector
v = norm(V); % magnitute of velcoity vector
vr = dot(R,V)/r; % radial component of velocity vector
if vr > 0
    fprintf('\n Satellite is moving away from perigee')
else
    fprintf('\n Satellite is moving towards the perigee')
end
H = cross(R,V); % h vector
h = norm(H); % magnitutde of h vector

inc = acosd(H(3)/h); % inclination angle

N = cross([0 0 1], H); % node line vector
n = norm(N); % magnitude of node line vector

if N(2) >= 0
    ra = acosd(N(1)/n); % right ascension
else
    ra = 360 - acosd(N(1)/n);
end

E = ((((v*v) - (mu/r))*R) - ((vr*r)*V))/mu; % eccentricity vector
e = norm(E); % magnitude of eccentricity vector

if E(3) >= 0
    w = acosd(dot(N, E)/(n*e));
else
    w = 360 - acosd(dot(N, E)/(n*e));
end

if vr >= 0
    ta = acosd(dot(E, R)/(e*r));
else
    ta = 360 - acosd(dot(E, R)/(e*r));
end
a = h*h/(mu*((e*e) - 1));
end