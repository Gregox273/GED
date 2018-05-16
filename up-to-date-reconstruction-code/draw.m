function draw(data, map);

% DRAW Draw an image
% 
%  DRAW(D, M) displays D as an image with colormap M.
%  If M is not given, greyscale is assumed.

% ensure we have the right input parameters
narginchk(1,2);
if (nargin==1)
  map = gray;
end

% draw image in current figure
hold off;

% images are usually defined the other way up
pcolor(flipud(data));

% change this to something else if you don't want linear interpolation
shading interp;

axis image;  % equal aspect ratio
axis off;       % no axes
set(gca,'Position',[0 0 1 1]); % fill the figure
set(gcf,'Color',[1 1 1]);        % white surround
colormap(map);
colorbar;
