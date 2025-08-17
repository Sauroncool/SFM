% known parameters
mu = 398600; 
delThorus = 2; % time (in hours) after which s/c state is needed
r0v = [20000, -105000, -19000]; % given initial position vector
v0v = [0.9, -3.4, -1.5]; % given initial velocity vector
r0 = norm(r0v); % r0 is magnitude of initial position
v0 = norm(v0v); % v0 is magnitude of initial velocity
tolerance = 10^-6;

% derived parameters
[e,h, thetha0] = orbital_parameters(r0v(1), r0v(2), r0v(3), v0v(1), v0v(2), v0v(3)); % calculating e, h and thetha0
vr0 = (dot(r0v, v0v))/r0; % magnitude of radial component of velocity vector
alpha = (2/r0) - (v0*v0/mu);
fprintf("\n")
if alpha < 0
    fprintf("Hyperbolic Trajectory")
elseif alpha > 0
    fprintf("Elliptical Trajectory")
else
    fprintf("Parabolic Trajectory")
end
a = 1/alpha; % semi major axis
delt = delThorus*3600; % time (in seconds) after which s/c state is needed

% Algorithm 3.3
x = sqrt(mu)*abs(alpha)*delt; % initial guess 

% coefficients of function f(x)
c0 = (r0*vr0)/sqrt(mu);
c1 = 1-(alpha*r0);
c2 = r0;
c3 = sqrt(mu)*delt;

ratio = 1; % initializing ratio
n = 1; % iterations
while abs(ratio) > tolerance
    z = alpha*x*x;
    f = c0*x*x*C(z) + c1*x*x*x*S(z) + c2*x - c3; % f
    fd = c0*x*((1-(alpha*x*x*S(z)))) + c1*x*x*C(z) + r0; % f' 
    ratio = f/fd;
    x = x - ratio;
    n = n + 1;
end

% lagrange coefficient and Algorithm 3.4
fL = 1 - (x*x*C(alpha*x*x)/r0);
gL = delt - (x*x*x*S(alpha*x*x)/sqrt(mu));
rv = fL*r0v + gL*v0v; % position vector of s/c after time t
r = norm(rv); % magnitude of position vector

fLd = sqrt(mu)*((alpha*x*x*x*S(alpha*x*x)) - x)/(r*r0);
gLd = 1 - (x*x*C(alpha*x*x)/r);
vv = fLd*r0v + gLd*v0v; % velocity vector of s/c after time t
v = norm(vv); % magnitude of velocity vector

% displaying outputs
fprintf("\nNumber of iterations= %d", n)
fprintf("\n")
fprintf("\nValue of universal variable= %d", x)
fprintf("\n")
fprintf("r = %d i + %d j + %d k", rv(1), rv(2), rv(3))
fprintf("\n")
fprintf("v = %d i + %d j + %d k", vv(1), vv(2), vv(3))
fprintf("\n")
