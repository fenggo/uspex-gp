function cart = Frac2Cart(frac, lat)
a     = lat(1);
b     = lat(2);
c     = lat(3);
alpha = lat(4);
beta  = lat(5);
gamma = lat(6);
v = sqrt(1 - cos(alpha)^2 - cos(beta)^2 -cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma));
M = [a b*cos(gamma) c*cos(beta); ...
0 b*sin(gamma) c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);...
0 0            c*v/sin(gamma)];
cart = (M*frac')';
