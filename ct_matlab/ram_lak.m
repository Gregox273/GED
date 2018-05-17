function [Y]=ram_lak(m, scale)

filter = -(m-1)/2:(m-1)/2;
w_max = pi/scale;

delta_w = 2*w_max/m;

filter = abs(filter*delta_w);

for i=1:m-1
    filter(i) = mean([filter(i), filter(i+1)]);
end

filter(m/2) = delta_w/4;

stairs(-(m-1)/2:(m-1)/2,filter);


Y = abs(filter)/(2*pi).*cos(filter./w_max*pi/2).^0.001;


stairs(-(m-1)/2:(m-1)/2, Y);