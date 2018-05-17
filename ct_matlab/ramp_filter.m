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

% apply filter to all angles
ram_lak_filter = ram_lak(m, scale);

ram_lak_filter = ones(angles,1)*ram_lak_filter;


%FFT input sinogram in the r direction for one angle. 

fast_fourier = fft(X,m,2); %If answer comes out wrong, change 2 to 1

fast_fourier = fftshift(fast_fourier,2);

filtered_signal = fast_fourier.*ram_lak_filter;

filtered_signal = fftshift(filtered_signal,2);

filtered = ifft(filtered_signal, [], 2);

Y = real(filtered(1:256, 1:256));

