function [V1, V2] = lamberts(R1, R2, t)
mu = 398600;
r1 = norm(R1);
r2 = norm(R2);
% prograde trajectory
R1_cross_R2 = cross(R1,R2);
if R1_cross_R2(3) >= 0
    deltaThetha = acosd(dot(R1,R2)/(r1*r2));
else
    deltaThetha = 360 - acosd(dot(R1,R2)/(r1*r2));
end

A = sqrt((r1*r2)/(1 - cosd(deltaThetha)))*sind(deltaThetha);

ratio = 1; % initializing ratio
n = 1; % iterations
tolerance = 10e-8;
z = 5;
while abs(ratio) > tolerance
    f = (((Y(r1, r2, z, A)/C(z))^1.5)*S(z)) + A*(Y(r1, r2, z, A))^0.5 - t*(mu)^0.5; % f
    if z == 0
        fd = sqrt(2)/40*Y(r1, r2, 0, A)^1.5 + A/8*(sqrt(Y(r1, r2, 0, A))+ A*sqrt(1/2/Y(r1, r2, 0, A)));
    else
        fd =  (Y(r1, r2, z, A)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z))+ 3*S(z)^2/4/C(z))+ A/8*(3*S(z)/C(z)*sqrt(Y(r1, r2, z, A))+ A*sqrt(C(z)/Y(r1, r2, z, A)));
    end
        
    ratio = f/fd;
    z = z - ratio;
    n = n + 1;
end

if z == 0
    fprintf('Parabolic Trajectory')
elseif z > 0
    fprintf('Elliptical Orbit')
else
    fprintf('Hyperbolic Trajectory')
end

lagf = 1 - Y(r1, r2, z, A)/r1;
lagg = A*(Y(r1, r2, z, A)/mu)^0.5;
laggd = 1 - Y(r1, r2, z, A)/r2;

V1 = (R2 - lagf*R1)/lagg;
V2 = (laggd*R2 - R1)/lagg;
end

