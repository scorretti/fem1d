function [x, y, varargout] = plotgrid (data)

% plotgrid                    Plot data on a grid
% 
% Description:                
%
% Input:
% - data                      Data provided as a list of points:
%                              - data(:,1) = x coordinate
%                              - data(:,2) = y coordinate
%                              - data(:,*) = any data to plot
%
% Output:                     
% - x, y                      Coordinates of the points on a grid, as provided by meshgrid()
% - ...                       Data organized on the grid
%
% Notes:                      
%
% Example:                    
% >> t = load('Ainstack.dat', 'ascii');
% >> [x, y, re, im] = plotgrid(t);
% >> figure ; contour(x, y, re, 'LineWidth', 2) ; axis image ; box on;
%
% See also:                   
%
% References:                 
%
% Validation:                 
%*
% Licence:                    Copyright Riccardo Scorretti
%                             This file is distributed under GPL-3.0-only ou GPL-3.0-or-later.
%
% Date:                       09-May-2017 - First version.

% --------------------------->| description of the function ---|------------------------------------------->| remarks

% Rebuild the structure of a grid
x_ = unique(data(:,1));       % = set of unique x coordinates
y_ = unique(data(:,2));       % = set of unique y coordinates

% Check that the data can be organized on a regular grid
nbPts = numel(x_)*numel(y_);
assert(nbPts == size(data,1), 'Data cannot be organized on a regular grid');

% Crate the grid
[x, y] = meshgrid(x_, y_);

nbData = size(data,2) - 2;    % = nb of quantities which can be plotted
for tx = 1 : numel(x_)
   for ty = 1 : numel(y_)
      % Retrieve the index (p) of the point at x_(tx), y_(ty)
      p = find((data(:,1) == x_(tx))  &   (data(:,2) == y_(ty)));
      for k = 1 : nbData
         varargout{k}(ty,tx) = data(p,2+k);
      end
   end
end

end


