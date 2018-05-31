function xtreme_reconstruct_all(H, FILENAME, ALPHA, material, METHOD)

% XTREME_RECONSTRUCT_ALL Create DICOM data from Xtreme RSQ file
%
%  XTREME_RECONSTRUCT_ALL(H, FILENAME, ALPHA) creates a series of DICOM
%  files for the Xtreme RSQ data given in the structure H. FILENAME is the
%  base file name for the data, and ALPHA is the power of the raised cosine
%  function used to filter the data.
%
%  XTREME_RECONSTRUCT_ALL(H, FILENAME, ALPHA, METHOD) can be used to
%  specify how the data is reconstructed. Possible options are:
%
%    'parallel' - reconstruct each slice separately using a fan to parallel
%                   conversion
%    'fdk' - approximate FDK algorithm for better reconstruction


% check inputs
narginchk(4,5);
if (nargin<5)
  METHOD = 'parallel';
end

% set frame number and DICOM UIDs for saving to multiple frames
z = 1;
seriesuid = dicomuid();
studyuid = dicomuid();
datetime = now();

% loop over each z-fan
for f=1:H.fan_scans:H.scans
  
  if (strcmp(METHOD,'fdk'))
    
    % correct reconstruction using FDK method would need to use all
    % h.fan_scans at once
    
  else
    
    % default method should reconstruct each slice separately
    for g=(f+H.skip_scans):(f+H.fan_scans-H.skip_scans-1)
      if (g<=H.scans)
        
        % reconstruct scan g:
        [f fmin fmax] = xtreme_get_rsq_slice(H, g);
        fmin = ones(278,1)*fmin;
        fmax = ones(278,1)*fmax;
        f = -log((f-fmin)./(fmax-fmin));  % Calibrate
        % Beam hardening correction:
        P = fake_source(material.mev,0.068,material.coeffs(:,15),6.7);
        % Get dimensions - air in ct_scan has depth 2*n*scale
        n = size(f, 2);
        % perform water calibration
        water = find(strcmp(material.name,{'Water'}));
        depth = [0:0.1:n/10-0.1];
        I_0_E = sum(P);
        water_scan = photons(P, material.coeffs(:, water), depth);
        I_tot = sum(water_scan);
        mu_tot= -log(I_tot/I_0_E);

        p1 = polyfit(mu_tot, depth, 4);
        estimated_depth = polyval(p1, mu_tot);

        p2 = polyfit(mu_tot(1:50), estimated_depth(1:50), 1);
        linear_depth = polyval(p2, mu_tot);

        % transform Y values
        for i =1:size(f,1)
            for j = 1:size(f,2)
                initial_mu = f(i,j);
                depth = polyval(p1, initial_mu);
                new_mu = depth/(p2(1)+10); % +10 is a fudge factor which seems to work well
                f(i,j) = new_mu;
            end
        end
        
        X = xtreme_fan_to_parallel(H, f);
        
        % Correct for zero column
        X = X(:,2:end);
        
        % Filter
        filtered_output = ramp_filter(X, H.scale/10, ALPHA);
        
        % Back project
        bp = back_project(filtered_output);

        % Hounsfield units
        %mu_w = 0.29;  % Empirical observation (of perspex tube)
        mu_w = 0.0785
        Y = (bp - mu_w)/mu_w * 1000;  % convert to HU

        % limit to normal DICOM range
        Y(Y<-1024) = -1024;
        Y(Y>3072) = 3072;
        
        draw(Y);
        % save the result as a DICOM file, with reference frame z
        create_dicom(Y, FILENAME, H.scale, H.scale, z, studyuid, seriesuid, datetime);

        % increment z
        z = z + 1;
        
      end
    end
    
  end
  
end