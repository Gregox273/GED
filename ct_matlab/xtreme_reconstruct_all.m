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
        X = xtreme_fan_to_parallel(H, f);
        % Correct for zero column
        X = X(:,2:end);
%         fmin = fmin(2:end);
%         fmax = fmax(2:end);
        
        % Filter
        filtered_output = ramp_filter(X, H.scale/10, ALPHA);
        
        % Back project
        bp = back_project(filtered_output);

        % Hounsfield units
        %Y = hu(, material, bp, H.scale/10);
        Y = bp;
        draw(Y);
        % save the result as a DICOM file, with reference frame z
        create_dicom(Y, FILENAME, H.scale, H.scale, z, studyuid, seriesuid, datetime);

        % increment z
        z = z + 1;
        
      end
    end
    
  end
  
end