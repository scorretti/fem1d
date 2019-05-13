function exemple_fem1d ()

% exemple_fem1d               
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

h = 0.025;                    % = mesh step

% Create the Fem1d object
x = linspace(0, 1, round(1/h));
FEM = Fem1d(x);
ide = 1 : FEM.getNbOfElements();

% Create the shape functions and their DOFs
U = ShapeFun1d('P1'); 
dofs_S = U.buildDofs(FEM, ide);

% Assembly a laplacian
coeff = @(x) 1 + 100*indfun(x, [0.5 1]);
WF = coeff*d(U)*d(U')

[A, b] = bilinear(coeff*d(U), d(U), FEM, ide, 'QUAD3');
[gA, gb] = FEM.assembly(A, b, dofs_S, dofs_S);

% Add Dirichlet boundary conditions
N = size(gA,1) ; gb = zeros(N, 1);
K = 1E6 * max(nonzeros(gA));                    % = huge constant used to impose Dirichlet contraints
gA(1,1) = gA(1,1) + K ; gb(1) = gb(1) + K*0;
gA(N,N) = gA(N,N) + K ; gb(N) = gb(N) + K*1;

% Solve the problem
u = gA \ gb;

% Plot the solution
figure
plot(x, u, 'o-');

end


