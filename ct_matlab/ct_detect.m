function [Y] = ct_detect(P, coeffs, depth, mas)

% CT_DETECT return detector photons for given material depths.
%
%  Y = CT_DETECT(P, coeffs, depth, mas) takes a source energy
%  distribution P (energies x 1), a set of material linear attenuation
%  coefficients coeffs (energies x materials), and a set of material depths
%  in depth (materials x samples) and returns the detections at each sample
%  in Y (1 x samples).
%
%  mas defines the current-time-product which affects the noise distribution
%  for the linear attenuation

% check inputs
narginchk(3, 4);
if (nargin<4)
  mas = 10000;
end

% extend source photon array so it covers all samples
Y = P*ones(1,size(depth,2));

% calculate array of residual mev x samples for each material in turn
for m=1:size(coeffs,2)
  Y = photons(Y, coeffs(:,m), depth(m,:), mas);
end

% sum over energy
Y = sum(Y);

% add in noise model
background_radiation = poissrnd(0.0001);
multiple_scattering_coefficient = 0.01;

no_source_photons = poissrnd(sum(P));

Y = Y + background_radiation + no_source_photons*multiple_scattering_coefficient;

% ensure it is above zero for log
Y(Y<=1) = 1;

