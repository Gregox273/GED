function [Y] = ct_calibrate(P, material, X, scale)

% CT_CALIBRATE convert CT detections to linearised attenuation
%
%  Y = CT_CALIBRATE(P, material, X, scale) takes the CT detection sinogram
%  in X (angles x samples) and returns a linear attenuation sinogram in Y
%  (angles x samples). P is the source energy distribution, material is the
%  material structure containing names, linear attenuation coefficients and
%  energies in mev, and scale is the size of each pixel in X, in cm.

% check inputs
narginchk(4,4);

% find coeffs corresponding to air
air = find(strcmp(material.name,{'Air'}));

% Get dimensions - air in ct_scan has depth 2*n*scale
n = size(X, 2);

% Perform calibration
air_scan = ct_detect(P, material.coeffs(:,air), 2*n*scale);

I_0_E = sum(air_scan);

Y= -log(X/I_0_E);


% perform water calibration
depth = [0:0.1:10];
I_0_E = sum(P);
water_scan = photons(P, material.coeffs(:, 5), depth);
I_tot = sum(water_scan);
mu_tot= -log(I_tot/I_0_E);

plot(mu_tot, depth)
xlabel('mu (cm^{-1})')
ylabel('depth (cm)')

%hold on

p1 = polyfit(mu_tot, depth, 4);
estimated_depth = polyval(p1, mu_tot);
%plot(mu_tot, estimated_depth)

p2 = polyfit(mu_tot(1:50), estimated_depth(1:50), 1);
linear_depth = polyval(p2, mu_tot);
%plot(mu_tot, linear_depth)


% transform Y values
for i =1:size(Y,1)
    for j = 1:size(Y,2)
        initial_mu = Y(i,j);
        depth = polyval(p1, initial_mu);
        new_mu = depth/p2(1);
        Y(i,j) = new_mu;
    end
end
