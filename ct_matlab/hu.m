function [Y] = hu(P, material, X, scale)

% HU convert CT reconstruction output to Hounsfield Units
%
%  Y = HU(P, material, X, scale) converts the output in X into Hounsfield
%  Units, using the material coefficients, photon energy P and scale given.

% check inputs
narginchk(4,4);

% find coeffs for water
water = find(strcmp(material.name,{'Water'}));

% convert to HU

% limit to normal DICOM range

