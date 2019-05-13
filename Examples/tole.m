function tole ()

% tole                        Solve the equation on a magnetic lamination
% 
% Description:                This program computes the magnetic field inside a magnetic lamination:
%
%                                  d      d
%                                  -- ( K -- T) + alpha (T - Ta) = 0
%                                  dx     dx
%
%                                T(0) = Te
%
%                                   d
%                                -K --T(L) = phi_e
%                                   dx
%
%                             This problem is written in weak form as:
%
%                                    d     d
%                                ( k --T , --T' ) + ( alpha T , T' ) - ( alpha Ta , T' ) + phi_e T'(L) = 0
%                                    dx    dx
%
%                             This program solves this problem, then compares the numerical solution with
%                             the analytical solution (1).
%
%
% Input:
% *                           
%
% Output:                     
% *                           
%
% Notes:                  1)   Matlab only, symbolic computation toolbox required.
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
%
% Date:                       12-Apr-2017 - First version.

% --------------------------->| description of the function ---|------------------------------------------->| remarks

close all ; clc

%% Parameters of the thermal model
m = 1 ; cm = 1E-2 ;mm = 1E-3; % = meter, centimeter, millimeter

th = 1*mm;                    % = thickness of the lamination

% Source term                 % = Magnetic field at the surface (A/m)
Hs = 1000;

% Physical parameters
sigma = 1E7;                  % = Electrical conductivity (S/m)
mur = 1000;                   % = Relative permeability (1)
fr = 1000;                      % = Frequency (hz)
omega = 2*pi*fr;              % = Pulsation (rad/s)


% Compute the skin depth
mu0 = pi*4E-7 ; mu = mu0*mur;
delta = sqrt(2/(omega*mu*sigma))


%% Create the FEM
N = 4;                       % = nb of elements
x = linspace(0, th/2, N);     % = coordinates of the nodes (m)

FEM = Fem1d(x); 
ide = FEM.getListOfElements();


%% Create the shape functions and associate it with the FEM
H = ShapeFun1d('P1');         % = temperature (C)
FEM.declare(H, ide);          % T is an unknown for the model FEM



%% Assembly the thermal problem

% To solve the problem without

% Assembly the FEM formulation
FEM.assembly(d(H)*d(H') + j*omega*sigma*mu*H*(H'), ide, 'QUAD4');
% FEM.assemblyBnd(phi_e*T', ide, 'right');     % <-- this is the term: phi_e*T'(L)

% Impose the constraint: T(0) = Te on the left boundary
FEM.impose(H == Hs, 'right');



%% Solve the problem
sol = FEM.solve();            % sol = values of the DOFs


% Plot the solution
figure
plot(x, abs(sol), 'o--') ; hold on
% H.plot(FEM, ide, linspace(0,1,10), sol, '-', 'MarkerSize', 6) ; hold on
% plot(x, Te-phi_e/K*x, 'r');


%% Plot the analytical solution
z = linspace(0, th/2, 1000);
a = (1+j)/delta;
H_ref = 1./sinh(a*th) * (Hs*sinh(a*(th/2+z)) + Hs*sinh(a*(th/2-z)));
plot(z, abs(H_ref));
legend('Finite Element', 'Analytical solution');
grid on
xlabel('z  (m)', 'FontSize', 16);
ylabel('|H|  (A/m)', 'FontSize', 16);

end


