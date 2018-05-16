function [Y] = ramp_filter(X, scale, alpha)

% RAMP_FILTER Ram-Lak filter with raised-cosine for CT reconstruction
%
%  Y = RAMP_FILTER(X) filters the input in X (angles x samples)
%  using a Ram-Lak filter.
%
%  Y = RAMP_FILTER(X, alpha) can be used to modify the Ram-Lak filter by a
%  cosine raised to the power given by alpha.

% check inputs
narginchk(2,3);
if (nargin<3)
  alpha = 0.001;
end

% get input dimensions
n = size(X,2);
angles = size(X,1);

% Set up filter length m to be a power of 2, and at least twice as long as input
% to prevent spatial domain aliasing and to speed up the FFT
m = floor(log(n)/log(2)+2);
m = 2^m;

% delete this %%%%
m = n;  % This makes it work, but what are we supposed to do with m?
%%%%%%%%%%%%%%%%%%

% apply filter to all angles
w_max = pi/scale;
delta_w = 2*w_max/(m-1);  % Spacing between frequency samples

% Matlab fft gives +ve frequencies in first half of vector, then -ve
% frequencies: need to rearrange filter to match this

filter = delta_w/4;  % Zero frequency value
w = 1:m/2;
not_zero_coeffs = w*delta_w;  % For positive frequencies
negative_coeffs = fliplr(not_zero_coeffs); % For negative frequencies
filter = cat(2, filter, not_zero_coeffs(1:end), negative_coeffs(2:end));

filter = abs(filter)/(2*pi) .* cos(filter/w_max * pi/2).^alpha;

% FFT of rows:
fast_fourier = fft(X,n,2);

% Turn filter into a matrix for element wise operation:
filter = ones(angles,1)*filter;
filter_output = fast_fourier.*filter;

Y = real(ifft(filter_output,n,2));

