classdef Jacobian1d < handle
   
   % Jacobian1d                  Transformation J : u=[0,1] --> x
   %
   % Description:                This class implements a quite simple jacobian, which takes the
   %                             reference domain [0,1] to Xn.
   %                             If no arguments are provided, the standard linear transformation is
   %                             used. However, the user can define its own transformation function:
   %
   %                                fx:(u,Xn) --> x
   %
   %                             where u are the coordinates with respect of the reference element [0 1],
   %                             and Xn=[Xn(1) Xn(2)] are the coordinates of the nodes of the real element.
   %
   %                             The computation of derivatives dx/du and du/dx is performed numerically.
   %
   % Properties:
   % - fx                        Handle to the function f(u,Xn) which computes the coordinates x in the
   %                             real space, given u = coordinates in the reference space, and Xn = 
   %                             coordinates of the nodes of the element in the real space.
   %
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
   % 10-Apr-2017 - First version.
   % 12-Apr-2017 - Improvement:  the methods x(), dxdu() and dudx() work for more than one element.
   % 14-Apr-2017 - Improvement:  the methods accepts Xn as well as (FEM, ide)
   
   % --------------------------->| description of the function ---|------------------------------------------->| remarks
   
   properties(GetAccess='public', SetAccess='protected')
      fx;                        % = handle to the function f : (u,Xn) --> x(u)
   end
   
   methods
      function J = Jacobian1d(fx)
         % Jacobian1d                  Constructors
         %
         % Description:
         %
         % Input: (copy constructor)
         % - J_                        Jacobian to be copied
         %
         % Output:
         % - J                         Deep copy of the Jacobian
         %
         %
         % Input:
         % - fx                        Handler to the function fx:(u,Xn) --> x 
         %                             (default = linear transformation)
         %
         % Output:
         % - J                         Jacobian
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 10-Apr-2017 - first version.
         
         if nargin == 1  &&  isa(varargin{1}, 'Jacobian1d')
            % Copy constructor
            J_ = varargin{1};    % = original object
            J.fx = J_.fx;
            return
         end
         
         if nargin == 0
            % Standard jacobian function = linear transformation which takes [0,1] to Xn
            J.fx = @(u,Xn) Xn(1)*(1-u) + Xn(2)*u;
         else
            % custom jacobian function
            J.fx = fx;
         end
      end
      
      function val = x(J, u, Xn)
         % x                           Compute the coordinates of x, given u
         %
         % Description:                This method computes the position x(u) for several quadrature points
         %                             and several elements.
         %                             The elements in the real space must be provided as a matrix Xn
         %                             of the form:
         %
         %                                Xn = [   x1(1)  x2(1)     <-- coordinates of the first element
         %                                         x1(2)  x2(2)
         %                                         .
         %                                         .
         %
         %                             The result is provided as a matrix val(e,q) where: e = index of the
         %                             element, q = index of the quadrature point.
         %
         % Input:
         % - J                         Jacobian
         % - u                         Coordinates in the reference domain (may be a vector)
         % - Xn                        Coordinates of the nodes of elements in the real space
         %
         % Output:
         % - val(e,q)                  Value of x(u)
         %
         % Notes:                  1)  The coordinates of the elements in the real space Xn can be provided
         %                             in two ways:
         %                              - Xn(e,2)     :  coordinates of the nodes
         %                              - {FEM, ide}  :  the method will retrieve the coordinates Xn(e,2)
         %                                               from the FEM object.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 10-Apr-2017 - first version.
         % 12-Apr-2017 - Improvement:  the methods x(), dxdu() and dudx() work for more than one element.
         % 14-Apr-2017 - Improvement:  the methods accepts Xn as well as (FEM, ide)
         
         if iscell(Xn)
            FEM = Xn{1} ; ide = Xn{2};
            [~ , Xn] = FEM.getElements('id', ide);
         end
         
         Ne = size(Xn,1);              % = nb of elements
         Nq = numel(u);                % = nb of quadrature points
         val = zeros(Ne, Nq);
         for e = 1 : Ne
            val(e,:) = feval(J.fx, u, Xn(e,:));
         end
      end
      function val = dxdu(J, u, Xn)
         % dxdu                        Compute the derivative dx/du, given u
         %
         % Description:                This method computes the derivative dx/du for several quadrature 
         %                             points and several elements.
         %                             The elements in the real space must be provided as a matrix Xn
         %                             of the form:
         %
         %                                Xn = [   x1(1)  x2(1)     <-- coordinates of the first element
         %                                         x1(2)  x2(2)
         %                                         .
         %                                         .
         %
         %                             The result is provided as a matrix val(e,q) where: e = index of the
         %                             element, q = index of the quadrature point.
         %
         % Input:
         % - J                         Jacobian
         % - u                         Coordinates in the reference domain (may be a vector)
         % - Xn                        Coordinates of the nodes of elements in the real space
         %
         % Output:
         % - val(e,q)                  Value of dx/du
         %
         % Notes:                  1)  The coordinates of the elements in the real space Xn can be provided
         %                             in two ways:
         %                              - Xn(e,2)     :  coordinates of the nodes
         %                              - {FEM, ide}  :  the method will retrieve the coordinates Xn(e,2)
         %                                               from the FEM object.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 10-Apr-2017 - first version.
         % 12-Apr-2017 - Improvement:  the methods x(), dxdu() and dudx() work for more than one element.
         % 14-Apr-2017 - Improvement:  the methods accepts Xn as well as (FEM, ide)
         
         if iscell(Xn)
            FEM = Xn{1} ; ide = Xn{2};
            [~ , Xn] = FEM.getElements('id', ide);
         end
         
         dfx = J.derive_u(J.fx);
         
         Ne = size(Xn,1);              % = nb of elements
         Nq = numel(u);                % = nb of quadrature points
         val = zeros(Ne, Nq);
         for e = 1 : Ne
            val(e,:) = feval(dfx, u, Xn(e,:));
         end
      end
      function val = dudx(J, u, Xn)
         % dudx                        Compute the derivative du/dx, given u
         %
         % Description:                This method computes the derivative du/dx for several quadrature 
         %                             points and several elements.
         %                             The elements in the real space must be provided as a matrix Xn
         %                             of the form:
         %
         %                                Xn = [   x1(1)  x2(1)     <-- coordinates of the first element
         %                                         x1(2)  x2(2)
         %                                         .
         %                                         .
         %
         %                             The result is provided as a matrix val(e,q) where: e = index of the
         %                             element, q = index of the quadrature point.
         %
         %                             The computation of du/dx is performed by using the inverse function
         %                             theorem, that is: given u such that dx/du != 0, one has:
         %
         %                                           1
         %                                du/dx =  -----
         %                                         dx/du
         %
         % Input:
         % - J                         Jacobian
         % - u                         Coordinates in the reference domain (may be a vector)
         % - Xn                        Coordinates of the nodes of elements in the real space
         %
         % Output:
         % - val(e,q)                  Value of du/dx
         %
         % Notes:                  1)  The coordinates of the elements in the real space Xn can be provided
         %                             in two ways:
         %                              - Xn(e,2)     :  coordinates of the nodes
         %                              - {FEM, ide}  :  the method will retrieve the coordinates Xn(e,2)
         %                                               from the FEM object.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 10-Apr-2017 - first version.
         % 12-Apr-2017 - Improvement:  the methods x(), dxdu() and dudx() work for more than one element.
         % 14-Apr-2017 - Improvement:  the methods accepts Xn as well as (FEM, ide)
         
         val = 1 ./ J.dxdu(u, Xn);
      end
   end
   
   methods(Static)
      function df = derive_u(f)
         % derive_u                    Given a function f(u,Xn), returns its (approximated) derivative df/du
         %
         % Description:                The derivative is approximated by using a centered finite difference,
         %                             with h = 1E-5.
         %                             If the function f(u,Xn) is NOT defined ouside the interval [0,1]
         %                             this function will fail utterly, and will have to be modified.
         %
         % Input:
         % - f                         Handle to a function f(u,Xn), with u in [0,1]
         %
         % Output:
         % - df                        Handle to a function df(u,Xn) which approximates df/du
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 10-Apr-2017 - first version.
         
         h = 1E-10;
         df = @(u,Xn) (f(u+h,Xn) - f(u-h,Xn))/(2*h);
      end
   end
end


