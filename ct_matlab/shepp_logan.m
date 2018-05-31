function [Y]=shepp_logan(m, scale)

filter = -m/2:(m-2)/2;
p_max = pi/scale;  % 1/(sample separation distance) = 2*f_max

delta_p = 2*p_max/m;  

filter = filter*delta_p

Y = abs(sin(filter*pi/(2*p_max)));

plot(filter, Y)
