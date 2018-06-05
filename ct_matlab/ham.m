function [Y]=ham(m, scale, alpha)
filter = -m/2:(m-2)/2;
w_max = pi/scale;  % 1/(sample separation distance) = 2*f_max

delta_w = 2*w_max/m;  % Spacing between discrete frequencies

filter = abs(filter*delta_w);  % Coefficients for w != 0

filter(m/2+1) = delta_w/4;  % DC value does not follow the above equation

Y = abs(filter)/(2*pi).*cos(filter/w_max*pi/2).^alpha;  % Raised cosine

window = hamming(m);

Y = Y.*window';