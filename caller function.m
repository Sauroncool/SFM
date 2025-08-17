%R1 = [5000 10000 2100];
%R2 = [-14600 2500 7000];
%t = 3600;
R1 = [5644 -2830 -4170];
R2 = [-2240 7320 -4980];
t = 20*60;
[V1, V2] = lamberts(R1, R2, t);
fprintf('\n')
fprintf('\n Solution of Lamberts Problem')
[h, e, ra, inc, w, ta, a] = sv_to_oe(R1, V1);
fprintf('\n Velocity vector at R1 (V1) = %d i + %d j + %d k', V1(1), V1(2), V1(3))
fprintf('\n Velocity vector at R2 (V2) = %d i + %d j + %d k', V2(1), V2(2), V2(3))
fprintf('\n')
fprintf('\n Orbital Elements')
fprintf('\n specific angular momentum (h) = %d km^2/s', h)
fprintf('\n eccentricity (e) = %d', e)
fprintf('\n Right Ascension of Ascending Node (RA) = %d degrees', ra)
fprintf('\n Argument of perigee (w) = %d degrees', w)
fprintf('\n True Anomaly (Thetha) = %d degrees', ta)
fprintf('\n Inclination (i) = %d degrees', inc)
fprintf('\n Semi major axis (a) = %d km', a)
fprintf('\n')