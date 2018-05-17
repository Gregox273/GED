function [Y]=ram_lak(m)

filter = -(m-1)/2:(m-1)/2;
w_max = pi;

delta_w = 2*w_max/m;

filter = abs(filter*delta_w);

for i=1:m-1
    filter(i) = mean([filter(i), filter(i+1)]);
end

filter = filter(1:m-1)

stairs(-(m-1)/2:(m-1)/2-1,filter);


Y = abs(filter)/(2*pi).*cos(filter./w_max*pi/2).^0.1;


stairs(-(m-1)/2:(m-1)/2-1, Y);