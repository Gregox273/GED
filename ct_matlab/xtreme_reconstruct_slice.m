function [Y] = xtreme_reconstruct_slice(FILENAME, slice)

h = xtreme_get_rsq_header(FILENAME)
[f fmin fmax] = xtreme_get_rsq_slice(h,slice);
p = xtreme_fan_to_parallel(h,f);
        
fmin = ones(278,1)*fmin;
fmax = ones(278,1)*fmax;
f = -log((f-fmin)./(fmax-fmin));  % Calibrate
X = xtreme_fan_to_parallel(h, f);
        
% Correct for zero column
X = X(:,2:end);
        
% Filter
filtered_output = ramp_filter(X, h.scale/10, 0.001);
        
% Back project
bp = back_project(filtered_output);

% Hounsfield units
mu_w = 0.29;  % Empirical observation (of perspex tube)
        
Y = (bp - mu_w)/mu_w * 1000;  % convert to HU

% limit to normal DICOM range
Y(Y<-1024) = -1024;
Y(Y>3072) = 3072;
        
draw(Y);
% save the result as a DICOM file, with reference frame z
%create_dicom(Y, h.filename, h.scale, h.scale, 100);
