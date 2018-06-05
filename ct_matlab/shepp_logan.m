function [Y]=shepp_logan(m, scale, alpha)
filter = -m/2:(m-2)/2;
w_max = pi/scale;  % 1/(sample separation distance) = 2*f_max

delta_w = 2*w_max/m;  % Spacing between discrete frequencies

filter_2 = abs(filter*delta_w);  % Coefficients for w != 0

filter_2(m/2+1) = delta_w/4;  % DC value does not follow the above equation

Y = abs(filter_2)/(2*pi).*cos(filter_2/w_max*pi/2).^alpha;  % Raised cosine

window = sinc(filter/m);

Y = Y.*window;