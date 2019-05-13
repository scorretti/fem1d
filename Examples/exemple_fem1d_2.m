function exemple_fem1d_2 ()

% exemple_fem1d_2               
% 
% Description:                
%
% Input:
% *                           
%
% Output:                     
% *                           
%
% Notes:                      
%
% Example:                    
%
% See also:                   
%
% References:                 
%
% Validation:                 
%
% Licence:                    Copyright Riccardo Scorretti
%                             This file is distributed under GPL-3.0-only ou GPL-3.0-or-later.
%
% Date:                       12-Apr-2017 - First version.

% --------------------------->| description of the function ---|------------------------------------------->| remarks

close all ; clc

h = 0.1;                    % = mesh step

% Define the properties of the "material"
coeff = @(x) 1 + 1*indfun(x, [0.5 1]);

% Create the FEM
x = linspace(0, 1, round(1/h));
FEM = Fem1d(x);
ide = 1 : FEM.getNbOfElements();

% Create the shape functions and associate it with the FEM
U = ShapeFun1d('P1'); 
FEM.declare(U, ide);

% Assembly a laplacian
% FEM.assembly(coeff*d(U)*d(U'), ide, 'QUAD3');

% FEM.assembly(coeff*d(U)*d(U'), 1:3, 'QUAD3');
% FEM.assembly(coeff*d(U)*d(U'), 4:numel(ide), 'QUAD3');


FEM.assembly(coeff*d(U)*d(U') + 1*(U'), ide, 'QUAD3');

% Impose the constraints
FEM.impose(U == 0, '-');
FEM.impose(U == 4, '+');

% Solve the problem
[u, lambda_] = FEM.solve();

% Plot the solution
figure
plot(x, u, 'o-');

end


