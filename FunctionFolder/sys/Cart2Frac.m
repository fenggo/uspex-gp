function frac = Cart2Frac(cart, lat)
a     = lat(1);
b     = lat(2);
c     = lat(3);
alpha = lat(4);
beta  = lat(5);
gamma = lat(6);
v = a*b*c*(1-cos(alpha)^2 - cos(beta)^2 -cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma))^(1/2);
M = [ 1/a, -cos(gamma)/(a*sin(gamma)), (b*cos(gamma)*c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)-b*c*cos(beta)*sin(gamma))/v;...
0,    1/(b*sin(gamma)),           -a*c*(cos(alpha)-cos(beta)*cos(gamma))/(v*sin(gamma)); ...
0,            0,                  a*b*sin(gamma)/v];
frac = (M*cart')';
