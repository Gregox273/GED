function [Y] = photons(X, coeff, depth, mas)

% PHOTONS calculate residual photons for a particular material and depth
%
%  Y = PHOTONS(X, coeff, depth, mas) takes the original energy in X (energy
%  x samples) and works out the residual energy Y (energy x samples) for a
%  particular material with linear attenuation coefficients given by coeff,
%  and a set of depths given by depth (1 x samples)
%
%  It is more efficient to calculate this for a range of samples rather then
%  one at a time
%
%  mas defines the current-time-product which affects the noise distribution
%  for the linear attenuation

% check inputs

% Work out residual energy for each depth and at each energy
mu_x = coeff*depth;

e_mu_x = exp(-(mu_x));

Y = X.*e_mu_x;

%Y = poissrnd(Y_one);

