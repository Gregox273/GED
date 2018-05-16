function xtreme_reconstruct_all(H, FILENAME, ALPHA, METHOD)

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
narginchk(3,4);
if (nargin<4)
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
        
        % reconstruct scan g

        % save the result as a DICOM file, with reference frame z

        % increment z
        z = z + 1;
        
      end
    end
    
  end
  
end