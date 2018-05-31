function [Y]=hamming(m, scale)

filter = -m/2:(m-2)/2;

p_max = pi/scale;  % 1/(sample separation distance) = 2*f_max

delta_p = 2*p_max/m; 

Y = 0.54-0.46*cos(2*pi*filter/(m-1));

plot(filter, Y)
