function [HU] = hu(P, material, X, scale)

% HU convert CT reconstruction output to Hounsfield Units
%
%  Y = HU(P, material, X, scale) converts the output in X into Hounsfield
%  Units, using the material coefficients, photon energy P and scale given.

% check inputs
narginchk(4,4);

% find coeffs for water
water = find(strcmp(material.name,{'Water'}));
% find coeffs for air
air = find(strcmp(material.name,{'Air'}));

n = size(X, 2);

water_scan = ct_detect(P, material.coeffs(:,water), 1); % per 1cm of water
air_scan = ct_detect(P, material.coeffs(:,air), 2*n*scale);

mu_w = -log(water_scan/air_scan);

% convert to HU
HU = (X - mu_w)/mu_w * 1000;

% limit to normal DICOM range
HU(HU<-1024) = -1024;
HU(HU>3072) = 3072;
