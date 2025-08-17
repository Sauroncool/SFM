G = 6.67e-20;  % km3kg-s-2
Re = 6378;  % km radius of earth
Rm = 1737;  % km radius of moon
rm = [384400 0 0];  % km vector joining earth to moon
D = norm(rm); % distance b/w earth and moon
mue = 398600.6698; % km3s-2
mum = 4902.8;  % km3s-2
vmvector = [0 1.0183 0]; % moon's velocity around earth
vm = norm(vmvector);
alpha0 = 28;
gamma0 = 6;
lamda = 55;
aepo = 320; % radius of EPO

% pg1
r0 =aepo + Re;
va = sqrt(mue/r0);
r0vector = -r0*[cosd(alpha0) sind(alpha0) 0];
ur0 = r0vector/r0;
Rs = 66183;
r2 = Rs;
r2vector = -Rs*[cosd(lamda) -sind(lamda) 0];
ur2 = r2vector/r2;
r1vector = rm + r2vector;
r1 = norm(r1vector);
ur1 = r1vector/r1;
deltathetha = acosd(dot(ur0,ur1));

%pg2
h1 = sqrt(mue*r0*(1 - cosd(deltathetha))/(r0/r1 + sind(deltathetha)*
tand(gamma0) - cosd(deltathetha)));
f = 1 - (mue*r1*(1 - cosd(deltathetha)))/(h1*h1);
g = r0*r1*sind(deltathetha)/h1;
gdot = 1 - (mue*r0*(1 - cosd(deltathetha)))/(h1*h1);
v0vector = (r1vector - f*r0vector)/g;
v0 = norm(v0vector);
vr0 = dot(v0vector, ur0);
v1vector = (gdot*r1vector - r0vector)/g;
v1 = norm(v1vector);
vr1 = dot(v1vector, ur1);
e1vector = ((v0*v0 - mue/r0)*r0vector - r0*vr0*v0vector)/mue;
e1 = norm(e1vector);
p1vector = e1vector/e1;
w1vector = cross(r0vector, v0vector)/h1;
q1vector = cross(w1vector, p1vector);

%pg3
a1 = h1*h1/(mue*(1-e1*e1));
T1 = 2*3.14*sqrt(a1*a1*a1/mue);
thetha0 = acosd(dot(p1vector, ur0));
t0 = (T1/(2*3.14))*((2*atand(sqrt((1-e1)/(1+e1))*tand(thetha0/2))) - 
(e1*sind(2*atand(sqrt((1-e1)/(1+e1))*tand(thetha0/2)))));
t0 = 130.37;
thetha1 = thetha0 + deltathetha;
t1 = (T1/(2*3.14))*((2*atand(sqrt((1-e1)/(1+e1))*tand(thetha1/2))) - 
(e1*sind(2*atand(sqrt((1-e1)/(1+e1))*tand(thetha1/2)))));
t1 = 239370;
deltat1 = (t1 - t0);
omegam = vm/D;
moonleadangle = omegam*deltat1;

%pg4
v2vector = v1vector - vmvector;
v2 = norm(v2vector);
vr2 = dot(v2vector, ur2);
h2vector = cross(r2vector, v2vector);
h2 = norm(h2vector);

%pg5
e2vector = cross(v2vector, h2vector)/mum - ur2;
e2 = norm(e2vector);
p2vector = e2vector/e2;
thetha2 = 360 - acosd(dot(p2vector, ur2));
t2 = (h2*h2*h2/(mum*mum*(e2*e2 - 1)^1.5))*((e2*sinh(2*atanh(sqrt((-1+e2)
/(1+e2))*tand(thetha2/2)))) - (2*atanh(sqrt
((-1+e2)/(1+e2))*tand(thetha2/2))));
deltat2 = 0 - t2;
deltattotal = deltat1 + deltat2;
rp2 = h2*h2/(mum*(1+e2));
zp2 = rp2 - Rm;
vp2 = sqrt((1+e2)*mum/rp2);

%pg6
deltav2 = sqrt(mum/rp2) - vp2;
deltav1 = sqrt(va*va + v0*v0 - 2*va*v0*cosd(gamma0));
deltavt = abs(deltav1) + abs(deltav2);

fprintf('\n')
fprintf('\n Position vector R0 = %d i + %d j + %d k', r0vector(1), 
r0vector(2), r0vector(3))
fprintf('\n Position vector R1 = %d i + %d j + %d k', r1vector(1), 
r1vector(2), r1vector(3))
fprintf('\n Position vector R2 = %d i + %d j + %d k', r2vector(1), 
r2vector(2), r2vector(3))
fprintf('\n Velocity vector V0 = %d i + %d j + %d k', v0vector(1),
v0vector(2), v0vector(3))
fprintf('\n Velocity vector V1 = %d i + %d j + %d k', v1vector(1), 
v1vector(2), v1vector(3))
fprintf('\n Velocity vector V2 = %d i + %d j + %d k', v2vector(1), 
v2vector(2), v2vector(3))
fprintf('\n Mission duration = %d days', deltattotal/3600/24)
fprintf('\n Eccentricity of hyperbolic lunar trajectory = %d', e2)
fprintf('\n delta-v from EPO to LTT = %d km/s', deltav1)
fprintf('\n delta-v for TLI = %d km/s', deltav2)
fprintf('\n total delta-v = %d km/s', deltavt)
fprintf('\n')