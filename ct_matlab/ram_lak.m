function [Y]=ram_lak(m, scale, alpha)
% Generates Ram-Lak filter (frequency domain) in fft shifted form where dc
% component is in the middle of the vector with negative frequencies to the
% left
factor = 32;

m_new=m/factor;

w_max = pi/(scale*factor);  % 1/(sample separation distance) = 2*f_max

filter = -m_new/2:(m_new-2)/2;

delta_w = 2*w_max/m_new;  % Spacing between discrete frequencies

filter = abs(filter*delta_w);  % Coefficients for w != 0

filter(m_new/2+1) = delta_w/4;  % DC value does not follow the above equation

Y_wrong_dim = abs(filter)/(2*pi).*cos(filter/w_max*pi/2).^alpha;  % Raised cosine

Y=[zeros(1, (m-m_new)/2), Y_wrong_dim, zeros(1, (m-m_new)/2)];

plot(-m/2:(m-2)/2,Y);
