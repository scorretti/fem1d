classdef Fem1d < handle
   
   % Fem1d                       A simple 1D Finite Element solver for didactical purposes
   %
   % Description:                This is a simple class to perform 1D FEM computations. It is not at all
   %                             efficient for large computations.
   %
   % Properties:                 The mesh is represented by a set of nodes (X). The n-th element is thus
   %                             composed of the nodes [ X(n) X(n+1) ]. The property idr stores the
   %                             ID of regions corresponding to each element. That is:
   %
   % - X                         Coordinates of the nodes
   % - idr                       ID-number of the k-th element = [X(k) X(k+1)]
   % - nbDofs                    Total number of DOFs of the declared shape functions
   % - Dofs                      Table of DOFs for the declared shape functions
   %
   %
   % Notes:                  1)  Naming convention for variables:
   %                             - h, k, n are used as index of loops
   %                             - Variables with the name beginning with "i_" mean: "Index of" ***
   %                             - idn, ide, idr = ID-number of Nodes, Elements, Regions
   %                             - Variables the name of which ends with "_" are temporary variables
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
   % 08-Apr-2017 - First version.
   % 13-Apr-2017 - Modification: a link between Fem1d and ShapeFun1d is created.
   % 16-Apr-2017 - Modification: property solution added.
   % 16-Apr-2017 - Modification: property J added.
   
   % --------------------------->| description of the function ---|---------------------------------------->| remarks
   
   properties(GetAccess='public', SetAccess='protected')
      X;                         % = coordinates of the nodes of the domain
      J;                         % = jacobian
      idr;                       % = ID of regions for each element
      nbDofs = 0;                % = total number of DOFs (Lagrange multipliers are NOT included)
      Dofs = {};                 % = table of DOFs for each declared shape function
      varNames = {};             % = name of the variables (experimental feature)
   end
   
   properties(GetAccess='public', SetAccess='public')
      A;                         % = global stiffness matrix
      b;                         % = global force vector
      solution = [];             % = numerical solution (if any)                       !!! not yet used !!!
   end
   
   methods
      % Constructors
      function F = Fem1d(varargin)
         % Fem1d                       Constructor
         %
         % Description:                The constructors creates a FEM. It is possible to specify the
         %                             coordinates of the nodes and (optionally) a jacobian which maps the
         %                             reference element [0, 1] to real elements.
         %
         %                             At present time, a single jacobian is defined for all elements.
         %                             This should be modified in future versions.
         %
         %
         %
         % Input: (copy constructor)
         % - F_                        FEM to be copied
         %
         % Output:
         % - F                         Deep copy of the FEM
         %
         %
         % Input:
         % - X                         Set of nodes.
         % - J                         Jacobian (default = linear map)
         %
         % Output:
         % - F                         FEM with the provided nodes. All elements are assigned to the
         %                             region 1.
         %
         %
         % Input:
         % *
         %
         % Output:
         % - F                         Empty FEM (use the method addRegion() to define the geometry)
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
         % 09-Apr-2017 - first version.
         % 16-Apr-2017 - Modification: property solution added.
         % 16-Apr-2017 - Modification: property J added.
         
         if nargin == 1  &&  isa(varargin{1}, 'Fem1d')
            % Copy constructor
            F_ = varargin{1};
            F.X = F_.X;
            F.J = F_.J;
            F.idr = F_.idr;
            F.nbDofs = F_.nbDofs;
            F.Dofs = F_.Dofs;
            F.A = F_.A;
            F.b = F_.b;
            F.solution = F_.solution;
            return
         end
         
         if nargin >= 1  &&  isnumeric(varargin{1})
            F.X = varargin{1} ; F.X = F.X(:);
            F.idr = ones(numel(F.X)-1, 1);
            
         elseif nargin == 0
            F.X = [];
            F.idr = [];
            
         else
            error('Unimplemented feature');
            % TODO
         end
         
         % Set the jacobian which maps the reference element [0, 1] to real elements [x1, x2]
         if nargin >= 2  &&  isa(varargin{2}, 'Jacobian1d')
            F.J = varargin{2};
         else
            F.J = Jacobian1d();
         end
      end
      
      % Methods which deal with elements and regions
      function ide = getListOfElements(F)
         % getListOfElements           Return the ID of all existing elements
         %
         % Description:
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - ide                       ID of all elements
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
         % 09-Apr-2017 - first version.
         
         ide = 1 : numel(F.idr);
      end
      function idr = getListOfRegions(F)
         % getListOfRegions            Return the ID of all existing regions
         %
         % Description:
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - idr                       ID of all regions
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
         % 09-Apr-2017 - first version.
         
         idr = unique(F.idr);
      end
      function nb = getNbOfRegions(F)
         % getNbOfRegions              Return the total number of regions
         %
         % Description:
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - nb                        Number of regions
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
         % 09-Apr-2017 - first version.
         
         nb = numel(F.getListOfRegions());
      end
      function nb = getNbOfElements(F)
         % getNbOfElements             Return the total number of elements
         %
         % Description:
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - nb                        Number of elements
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
         % 09-Apr-2017 - first version.
         
         nb = numel(F.idr);
      end
      function [ide, X] = getElements(F, argName, argValue)
         % getElements                 Return the ID and the elements in the form X(:,2)
         %
         % Description:                This method search a subset of elements, and returns their ID and
         %                             a copy of the elements in the form X(:,2). That is, the n-th element
         %                             is composed of the nodes of coordinates [ X(n,1) X(n,2) ].
         %
         % Input:
         % - F                         FEM
         % - argName                   Criterion of search, which may be one of: 'id', 'region'
         % - argValue                  Argument of the criterion
         %
         % Output:
         % - ide                       ID of elements
         % - X                         Copy of the elements in the form: X(:,2)
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
         % 09-Apr-2017 - first version.
         
         % If required, search the elements
         switch lower(argName)
            case {'id' , 'ide'}
               ide = argValue;
            case 'region'
               ide = find(ismember(F.idr, argValue));
            otherwise
               error('Unimplemented criterion for searching elements');
         end
         
         % If the elements are not required, exit
         if nargout < 2 , return ; end
         
         % Generate an easy-to-use copy of the required elements
         X = zeros(numel(ide), 2);
         X(:,1) = F.X(ide);
         X(:,2) = F.X(ide+1);
      end
      function [idr, ide] = addRegion(F, X, idr)
         % addRegion                   Add a region
         %
         % Description:                This method allows to add a group of (contiguous) nodes, which
         %                             will belong to a given region. There is no constraint on the
         %                             arguments: new regions can intersect existing regions (in which case
         %                             elements will be assigned to the new region).
         %
         % Input:
         % - F                         FEM
         % - X                         Nodes to add to the FEM
         % - idr                       ID of the region (can be assigned by default)
         %
         % Output:
         % - idr                       ID of the region
         % - ide                       ID of the elements which compose the region
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
         % 09-Apr-2017 - first version.
         
         % Ensure that the nodes are sorted
         X = sort(X) ; X = X(:);
         
         % If no ID-number is provided for the region, assign to it a new ID-number by default
         if nargin < 3 , idr = F.getNextFreeID(F.idr) ; end
         
         % Merge the existing set of nodes with the existing nodes
         X_  = union(F.X, X);
         idr_ = zeros(numel(X)-1, 1);
         
         % Assign to each each element the right region, which is determined by its first node.
         [~, t_] = ismember(F.X(1:end-1), X_(1:end-1));
         idr_(t_) = F.idr;
         [~, t_] = ismember(X(1:end-1), X_(1:end-1));
         idr_(t_) = idr;
         ide = t_;
         
         F.X = X_(:);
         F.idr = idr_(:);
         
         F.check();
      end
      function J = getJacobian(F, ide)
         % getJacobian                 Return the jacobian for the given elements
         %
         % Description:                This method returns the jacobian which maps the reference element
         %                             [0, 1] to real elements.
         %
         %                             All the elements must have the same jacobian, otherwise an error
         %                             will be emitted. At present time, this constraint is satisfied
         %                             by construction (this will be modified in future versions).
         %
         % Input:
         % - F                         FEM
         % - ide                       ID of the elements
         %
         % Output:
         % - J                         Jacobian which maps [0, 1] to real elements
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
         % 09-Apr-2017 - first version.
         J = F.J;
      end
      
      % Methods which deal with variables and DOFs
      function id = declare(F, S, ide)
         % declare                     Declare as unknown a variable defined by its shape function
         %
         % Description:                This method creates a link between a shape function and the FEM.
         %                             More precisely, this method:
         %                              - assign a new ID-number to the shape function, and creates the
         %                                link FEM <--> shape function,
         %                              - builds the table of DOFs numbering for the shape function.
         %
         %                             After a variable has been declared by a FEM, it is possible to
         %                             retrieve the corresponding DOFs numbering.
         %
         %                             A single variable can be associated with many FEMs. The DOFs
         %                             numbering is stored on the side of the FEM. This is consistent with
         %                             the fact that the number of DOFs depends on the elements over which
         %                             a variable is defined: therefore, DOFs numbering is an intrinsic
         %                             property of the mesh (i.e. of the FEM), not of the shape function.
         %
         %                             The numberig of the DOFs of the same variable will generally be
         %                             different from a FEM object to another.
         %
         %                             If possible, the name of the variable S is retrieved and stored in
         %                             the property varName.
         %
         % Input:
         % - F                         FEM
         % - S                         Shape function
         % - ide                       ID-number of elements where the shape function has a value
         %                             (default = all elements)
         %
         % Output:
         % - id                        ID-number of the variable, with respect of this FEM
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
         % 13-Apr-2017 - first version.
         % 19-Apr-2017 - Modification: property varName added (experimental).
         
         % By default, shape functions are defined over all of the elements
         if nargin < 3 , ide = 1 : F.getNbOfElements() ; end
         
         % Assign an ID-number to the shape function and associate it with the FEM
         id = numel(F.Dofs) + 1;
         S.associate(F, id);
         
         % Generate the table of DOFs (dofs_) for the selected elements (ide)
         [dofs_, nb] = S.buildDofs(F, ide);     % = table of DOFs, total nb of DOFs generated
         
         % Create a table of DOFs (dofs) which accounts all elements. Copy the entries of dofs_ to the
         % corresponding entries of dofs. The lines corresponding to elements which take no DOF are
         % filled with 0 (this could be optimized to save memory).
         % An offset is added to each DOF number, so that all DOFs have a different number, and DOF
         % numbering is consecutive
         dofs = zeros(F.getNbOfElements, size(dofs_,2));
         dofs(ide,:) = dofs_ + F.nbDofs;
         
         % Save the table of DOFs with the ID (id) corresponding to the shape function (S) and increase the
         % total number of DOFs
         F.Dofs{id} = dofs;
         F.nbDofs = F.nbDofs + nb;
         
         % Try to save the name of the variables
         F.varNames{id} = inputname(2);
      end
      function dofs = getDofs(F, id)
         % getDofs                     Return the table of DOFs for a given variable
         %
         % Description:                This method returns the table of DOFs for a given variable, or for
         %                             a variable with a given ID-number. The table of DOFs is a matrix
         %                             of the form dofs(e,p), where e = index of elements, p = index of DOFs.
         %
         %                             The table of DOFs is defined over all of the elements. Elements where
         %                             the variable is not defined take the DOF ID-number 0, meaning that
         %                             they are not associated with that quantity.
         %
         % Input:
         % - F                         FEM
         % - id                        ID-number of the variable
         %
         %
         % Input:
         % - F                         FEM
         % - S                         Variable (i.e. shape function)
         %
         % Output:
         % - dofs(e,p)                 Table of the DOFs numbers over all of the elements
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
         % 13-Apr-2017 - first version.
         
         % If a shape function is provided, retrieve its ID-number (for this particular FEM)
         if isa(id, 'ShapeFun1d') , id = S.getId(F) ; end
         
         dofs = F.Dofs{id};
      end
      function [A, b] = getLinearSystem(F)
         % getLinearSystem             Get the linear system Ax = b
         %
         % Description:                This method returns the global linear system of the FEM. If the
         %                             linear system has not yet been read (that is, if A = [] and b = [])
         %                             then this method creates an empty matrix A and vector b.
         %
         %                             This method must be called AFTER all variables have been declared
         %                             by using the method declare().
         %
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - A                         Global stiffness matrix
         % - b                         Global force vector
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
         % 13-Apr-2017 - first version.
         
         if isempty(F.A)
            F.A = sparse(F.nbDofs, F.nbDofs);
            F.b = zeros(F.nbDofs, 1);
         end
         A = F.A ; b = F.b;
      end
      function reset(F, opt)
         % reset                       Reset the FEM
         %
         % Description:
         %
         % Input:
         % - F                         FEM
         % - opt                       Possible options:
         %                              - 'all'    = clear the linear system and the DOFs (default)
         %                              - 'system' =   "    "    "     "     only
         %
         % Output:
         % *
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
         % 13-Apr-2017 - first version.
         
         if nargin < 2 , opt = 'all' ; end
         switch lower(opt)
            case 'all'                 % Clear the linear system and the DOFs
               F.nbDofs = 0;
               F.Dofs = {};
               F.A = [];
               F.b = [];
               F.solution = [];
               
            case 'system'              % Clear the linear system only
               F.A = [];
               F.b = [];
               F.solution = [];
               
            otherwise
               error('Unimplemented feature');
         end
      end
      function [x, lambda] = solve(F)
         % solve                       Solve the FEM problem
         %
         % Description:                This method solves the linear system of the FEM, stores the
         %                             full solution (aka, together with Lagrange multipliers, if any) and
         %                             returns separately the value of the DOFs and of Lagrange multipliers.
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % - x                         Value of the DOFs
         % - lambda                      "   "  Lagrange multipliers (if any)
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
         % 13-Apr-2017 - first version.
         
         % Solve the linear system
         u = F.A \ F.b;
         x = u(1:F.nbDofs);
         lambda = u(F.nbDofs+1:end);
         
         % Store the solution
         F.solution = u;
      end
      function hndl = spy(F, varargin)
         % spy                         Display the sparsity pattern of the stiffness matrxi
         %
         % Description:                This method is an improved version of the standard function spy(),
         %                             in that it display the boundaries between the different submatrix
         %                             which compose the global stiffness matrix, and displays (if possible)
         %                             the names of the unknowns of the problem.
         %
         % Input:
         % - F                         FEM
         % - varargin                  Additional arguments to pass to the function spy()
         %
         % Output:
         % - hndl                      Handles to the lines created
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
         % 19-Apr-2017 - first version.
         
         % First, invoke the standard spy() function
         spy(F.A, varargin{:}) ; hold on
         
         % Draw the boundaries of the submatrix of A
         N = size(F.A, 2) ; hndl = [];
         xt = [] ; xtlab = {};
         for k = 1 : numel(F.Dofs)
            t1 = min(setdiff(F.Dofs{k}(:), [0]));
            t2 = max(F.Dofs{k}(:));
            hndl(end+1) = plot([0 N], [t2 t2]+0.5, 'k--');
            hndl(end+1) = plot([t2 t2]+0.5, [0 N], 'k--');
            
            % Store the values for X and Y ticks
            xt(end+1) = t1 ; xtlab{end+1} = num2str(t1);
            xt(end+1) = (t1+t2)/2 ; xtlab{end+1} = F.varNames{k};
         end
         axis image
         
         if t2 < N
            xt(end+1) = t2+1 ; xtlab{end+1} = num2str(t2+1);
            xt(end+1) = (t2+1+N)/2 ; xtlab{end+1} = '\lambda';
         end
         xt(end+1) = N ; xtlab{end+1} = num2str(N);
         
         % Set the ticks to display the name of the variables
         set(gca, 'XTick', xt(2:2:end));
         set(gca, 'XTickLabel', xtlab(2:2:end));
         set(gca, 'YTick', xt(1:2:end));
         set(gca, 'YTickLabel', xtlab(1:2:end));
         set(gca, 'XAxisLocation', 'top');
      end
      
      % Assembly routines
      function [gA, gb] = assembly(F, varargin)
         % assembly                    Assembly equations to the global linear system
         %
         % Description:                This method dispatch the execution to the right subroutine for
         %                             assemblying equations. There are several possibilities:
         %                              1) assembly "raw" linear equations, provided as element-wise
         %                                 matrix A(e,1,i,j) and b(e,1,1,j).
         %                              2)  assembly one or more bilinear forms
         %                              3)  assembly constraints               !!! unimplemented feature !!!
         %
         %                             In fact, in the case 1) the argument FEM is not used. The user
         %                             must provide the global system gA.x = gb of the appropriate size.
         %                             This subroutine is used by the other subroutine, more user-friendly.
         %
         %
         % Input (assembly raw equations):
         % - F                         FEM
         % - A(e,1,i,j)                Element-wise stiffness matrix
         % - b(e,1,1,j)                   "     "   force vector
         % - dof_i(e,i)                DOFs numbering for the unknown       (indexed by i)
         % - dof_j(e,j)                 "      "       "   "  test function (   "    "  j)
         % - gA                        Global stiffness matrix
         % - gb                          "    force vector
         %
         % Input (assembly bilinear forms):
         % - F                         FEM
         % - B                         Bilinear form(s)
         % - ide                       ID-number of the elements where bilinear forms must be evaluated
         % - u                         Coordinates of the computation points, or quadrature formula
         % - flag                      If true (default), update the linear system of the FEM
         %
         % Output:
         % - gA                        Global stiffness matrix
         % - gb                        Global force vector
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
         % 13-Apr-2017 - first version.
         
         if isnumeric(varargin{1})
            [gA, gb] = F.assembly_1(varargin{:});
            
         elseif isa(varargin{1}, 'Bilinear1d')
            [gA, gb] = F.assembly_2(varargin{:});
            
         else
            error('Unimplemnted feature (this may be due to a missing ('') to mark the test function)');
         end
      end
      function [gA, gb] = assembly_2(F, B, ide, u, flag)
         % assembly_2                  Assembly bilinear forms to the global linear system
         %
         % Description:                This method evaluates bilinear forms (B) in a subset of elements
         %                             (ide). The local coordinates may be provided directly (u = u(q)),
         %                             or through a quadrature formula (u = 'QUAD3').
         %
         %                             The linear system which is considered is the one associated with
         %                             the FEM. However, the user has the possibility of not updating
         %                             the linear system of FEM.
         %
         % Input:
         % - F                         FEM
         % - B                         Bilinear form(s)
         % - ide                       ID-number of the elements where bilinear forms must be evaluated
         % - u                         Coordinates of the computation points, or quadrature formula
         % - flag                      If true (default), update the linear system of the FEM
         %
         % Output:
         % - gA                        Global stiffness matrix
         % - gb                        Global force vector
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
         % 13-Apr-2017 - first version.
         
         % By default, use the linear transformation
         J = F.getJacobian(ide);
         
         % Bay default, update the linear system of the FEM
         if nargin < 5 , flag = true ; end
         
         [gA, gb] = F.getLinearSystem();
         for k = 1 : numel(B)
            [dof_i, dof_j] = B(k).getDofs(F);   % = DOFs associated with the lhs, rhs of the bilinear form
            [A, b] = B(k).feval(F, ide, u);     % = element-wise matrix corresponding to the bilinear form
            
            if isempty(dof_i)
               [gA, gb] = F.assembly_1(A, b,  [], dof_j(ide,:), gA, gb);
            else
               [gA, gb] = F.assembly_1(A, b, dof_i(ide,:), dof_j(ide,:), gA, gb);
            end
         end
         
         % Update the linear system
         if flag
            F.A = gA ; F.b = gb;
         end
      end
      function [gA, gb] = assembly_1(F, A, b, dof_i, dof_j, gA, gb)
         % assembly                    Assembly routine
         %
         % Description:                This method assemblies element-wise stiffness matrix A(e,1,i,j) and
         %                             force vector b(e,1,1,j) into a global linear system gA.x = gb
         %
         % Input (generic assembly routine):
         % - F                         FEM (this argument is actually not used)
         % - A(e,1,i,j)                Element-wise stiffness matrix
         % - b(e,1,1,j)                   "     "   force vector
         % - dof_i(e,i)                DOFs numbering for the unknown       (indexed by i)
         % - dof_j(e,j)                 "      "       "   "  test function (   "    "  j)
         % - gA                        Global stiffness matrix
         % - gb                          "    force vector
         %
         % Output:
         % - gA                        Updated global stiffness matrix
         % - gb                           "      "    force vector
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
         % 13-Apr-2017 - first version.
         
         if nargin < 6 , gA = [] ; gb = [] ; end
         
         % if A is not empty, assembly it into a sparse matrix
         if ~isempty(A)
            [Ne, ~, Ni, Nj] = size(A);             % = nb of elements (Ne), size of el. st. mat (Ni,Nj)
            i_ = zeros(Ne, Ni, Nj);
            j_ = zeros(Ne, Ni, Nj);
            
            % Fill the vectors i_ and j_ with the right number of DOFs
            for i = 1 : Ni
               for j = 1 : Nj
                  i_(:,i,j) = dof_i(:,i);
                  j_(:,i,j) = dof_j(:,j);
               end
            end
            
            % Remove the terms which must not be assembled
            t_ = find((i_ ~= 0)  &  (j_ ~= 0));
            i_ = i_(t_) ; j_ = j_(t_) ; A = A(t_);
            
            if isempty(gA)
               % No global stiffness matrix is provided
               gA = sparse(j_(:), i_(:), A(:));
            else
               % A global stiffness matrix is provided: accumulate the result into it
               gA = gA + sparse(j_(:), i_(:), A(:), size(gA,1), size(gA,2));
            end
         end
         
         % If b is not empty, assembly it
         if ~isempty(b)
            [Ne, ~, ~, Nj] = size(b);              % = nb of elements (Ne), size of force vec. (1,Nj)
            i_ = ones(Ne, 1, Nj);
            j_ = zeros(Ne, 1, Nj);
            
            % Fill the vector j_ with the right number of DOFs
            for j = 1 : Nj
               j_(:,1,j) = dof_j(:,j);
            end
            
            % Remove the terms which must not be assembled
            t_ = find(j_ ~= 0);
            i_ = i_(t_) ; j_ = j_(t_) ; b = b(t_);
            
            % REMARK: the equation is computed in the form: Ax+b = 0, therefore the linear system
            % writes Ax = -b
            if isempty(gb)
               % No global force vector is provided
               gb = -full(sparse(j_(:), i_(:), b(:)));
            else
               % A global force vector is provided: accumulate the result into it
               gb = gb - full(sparse(j_(:), i_(:), b(:), size(gb,1), 1));
            end
         end
      end
      function [gA, gb] = assemblyBnd(F, B, ide, bnd, varargin)
         % assemblyBnd                 Assembly bilinear forms on a boundary
         %
         % Description:                This method evaluates bilinear forms (B) on the boundary of a
         %                             subset of elements, which are assumed to be contiguous.
         %
         %                             The possible options are:
         %                              - '+', or 'right'      : compute the bilinear form on the right boundary
         %                              - '-', or 'left'       :    "     "     "      "   "   "  left    "
         %                              - ']', or 'jump-right' : compute the jump of the bilinear form
         %                                                       through the right boundary
         %                              - '[', or 'jump-left'  : compute the jump of the bilinear form
         %                                                       through the left boundary
         %
         %
         % Input:
         % - F                         FEM
         % - B                         Bilinear form(s)
         % - ide                       ID-number of the elements where bilinear forms must be evaluated
         % - nbd                       Boundary ('+', '-', ']', '[')
         % ...                         Additional arguments passed to assembly()
         %
         % Output:
         % - gA                        Global stiffness matrix
         % - gb                        Global force vector
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
         % 16-Apr-2017 - first version.
         
         switch lower(bnd)
            case {'-', 'left'}         % left boundary
               [gA, gb] = assembly_2(F, -B, min(ide), 0, varargin{:});
               
            case {'+', 'right'}        % right boundary
               [gA, gb] = assembly_2(F, B, max(ide), 1, varargin{:});
               
            case {'[', 'jump-right'}   % jump on the right boundary
               [gA, gb] = assembly_2(F, B, max(ide), 1, varargin{:});
               [gA, gb] = assembly_2(F, -B, max(ide)+1, 1, varargin{:});
               
            case {'[', 'jump-left'}   % jump on the left boundary
               [gA, gb] = assembly_2(F, -B, min(ide), 1, varargin{:});
               [gA, gb] = assembly_2(F, B, min(ide)-1, 1, varargin{:});
               
            otherwise
               error('Boundary type (%s) unknown', bnd);
         end
      end
      
      % Methods for imposing constraints
      function [gA, gb] = impose(F, C, ide, u)
         % impose                      Modify the global linear system so as to impose constraints
         %
         % Description:                This method imposes the constraint (C) by using Lagrange
         %                             multipliers.
         %                             First, the bilinear form C is evaluated, so as to obtain a set of
         %                             linear equations:
         %
         %                                B.x = H
         %
         %                             Then the global linear system A.x = b is modified so as to impose
         %                             the constraint B.x = H by using the method of Lagrange multipliers,
         %                             that is:
         %
         %                                A.x = b        ---->       | A  B'|.| x      | = | b |
         %                                                           | B  0 | | lambda |   | H |
         %
         %
         % A tehcnical point:          In order to avoid programming again operators, constraints are
         %                             stocked as bilinear form with a "dummy" test function. For instance,
         %                             the constraint:
         %
         %                                A == 0
         %
         %                             will be represented by dummy bilinear form:
         %
         %                                (A, *)
         %
         %                             Therefore, C must be an instance of the class Bilinear1d.
         %
         %
         % Input:
         % - F                         FEM
         % - C                         Constraint
         % - ide                       ID-number of the elements where bilinear forms must be evaluated
         % - u                         Boundary where the constraint has to be imposed:
         %                              - '-' or 'left'  : left boundary
         %                              - '+' or 'right' : right   "
         %
         % Output:
         % - gA                        Global stiffness matrix
         % - gb                        Global force vector
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
         % 16-Apr-2017 - first version.
         
         assert(isa(C, 'Bilinear1d'));
         
         % This is a shortcut to simplify the synopsys
         if nargin < 4
            if ischar(ide)
               % Synopsys: impose(FEM, C, 'left') or impose(FEM, C, 'right')
               switch(ide)
                  case {'+', 'right'}
                     ide = F.getNbOfElements();
                     u = 1;
                  case {'-', 'left'}
                     ide = 1;
                     u = 0;
                  otherwise
                     error('Unimplemented feature');
               end
               
            elseif iscell(ide)
               % Synopsys: impose(FEM, C, {ide_1, ide_2})
               if isnumeric(ide{1})  &&  isnumeric(ide{2})
                  if max(ide{1}) == min(ide{2})-1
                     ide = [ max(ide{1}) min(ide{2}) ];
                     u = [1 0];
                  elseif max(ide{2}) == min(ide{1})-1
                     ide = [ max(ide{2}) min(ide{1}) ];
                     u = [1 0];
                  else
                     warning('A constraint could not be imposed');
                     return
                  end
               end
               
            else
               error('Unimplemented feature');
            end
         end
         
         assert(numel(ide) == numel(u));
         
         % Obtain the jacobian
         J = F.getJacobian(ide);
         
         % Evaluate the bilinear form to obtain the constraint equation Bx = H
         B = sparse(1,size(F.A,2)) ; H = [0];
         
         for e = 1 : numel(ide)
            for k = 1 : numel(C)
               [dof_i, ~] = C(k).getDofs(F);   % = DOFs associated with the lhs of the bilinear form
               [A_, b_] = C(k).feval(F, ide(e), u(e));
               if isempty(dof_i)
                  [B, H] = F.assembly_1(A_, b_, [], 1, B, H);
               else
                  [B, H] = F.assembly_1(A_, b_, dof_i(ide(e),:), 1, B, H);
               end
            end
         end
         
         % Add the equation B.x = H to the global linear system
         [F.A, F.b] = F.lagrange(B, H, F.A, F.b);
         [gA, gb] = F.getLinearSystem();
      end
      function [A, b] = lagrange(F, B, H, A, b, opt_Dirichlet)
         % lagrange                    Modify the linear system A.x = b so as to impose B.x = H
         %
         % Description:                This methods modifies the linear system A.x = b so as to impose
         %                             the constraint B.x = H by using Lagrange multipliers, that is:
         %
         %                                A.x = b        ---->       | A  B'|.| x      | = | b |
         %                                                           | B  0 | | lambda |   | H |
         %
         %
         % Input:
         % - F                         FEM
         % - B, H                      Equations of the constraint B.x = H
         % - A, b                      Linear system A.x = b
         % - opt_Dirichlet             Precises how to impose Dirichlet constraints:
         %                              - 0 = Lagrange multipliers
         %                              - 1 = penalty
         %                              - 2 = elimiation
         %
         % Output:
         % - A, b                      Modified linear system A.x = b
         %
         % Notes:                  1)  This function is for internal usage only.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 09-Apr-2017 - first version.
         % 19-Apr-2017 - Improvement:  it is checked if absurd constraints (i.e. 0 == *) are imposed.
         
         if nargin < 5 , A = F.A ; b = F.b ; end
         if nargin < 6 , opt_Dirichlet = 0 ; end
         
         [Nj, Ni] = size(B) ; assert(size(B,1) == Nj);
         
         % Remove empty lines from B and H
         t = [];
         for k = 1 : Nj
            if nnz(B(k,:)) > 0
               t(end+1) = k;
            elseif H(k) ~= 0
               warning('Constraint of the form 0 == * found');
            end
         end
         B = B(t,:) ; H = H(t);
         
         % Handle Dirichlet constraints
         for k = 1 : Nj
            if nnz(B(k,:)) == 1
               % The constraint to impose is: x(j) = vb/va
               [i, j, va] = find(B(k,:));
               vb = H(i);
               xj = vb/va;
               
               switch opt_Dirichlet                                                                         % TO DO !!!
                  case 0      % Lagrange multipliers (nothing to do)
                     ;
                  
                  case 1      % Penalty method
                     warning('Penalty method is not yet implemented');
                     ;
                  
                  case 2      % Elimination
                     % Modify b so as to enforce Dirichlet constraint
                     % 
                     [NAj, NAi] = size(A);
                     b = b - xj*A(:,j);
                     
                     % Modify the linear system A.x = b so as to enforce Dirichlet constraint
                     % A(:,j) = 0 ; A(j,:) = 0 ; A(j,j) = 1 ; b(j) = xj;
                     
                     % Eliminate the corresponding DOF and renumber the unknowns
                     % [A, b] = F.removeDof(j);
                     
                     % Eliminate the corresponding equation from B.x = H
                     ti_ = setdiff(1:Ni, i);
                     tj_ = setdiff(1:Nj, j);
                     B = B(tj_,ti_) ; H = H(tj_);
                     
                  otherwise
                     error('Unimplemented feature');
               end
            end
         end
         
         A = [ A B' ; B sparse(Nj,Nj) ];
         b = [ b ; H ];
         
         if nargout == 0 , F.A = A ; F.b = b ; end
      end
      
      % This method executes some consistency tests on the mesh
      function check(F)
         % check                       Perform some tests on the FEM (for internal usage only)
         %
         % Description:                This method performs some (in principle useless) tests the following
         %                             consistency condition:
         %                              - the number of nodes and elements are consistent
         %                              - the numbering of DOFs is correct
         %
         % Input:
         % - F                         FEM
         %
         % Output:
         % *
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
         % 09-Apr-2017 - first version.
         % 13-Apr-2017 - Modification: the numbering of DOFs is also checked.
         
         % Check that theare are no errors in the mesh
         assert( ...
            (isempty(F.X)  &&  isempty(F.idr))  ||         ...
            (numel(F.X) == numel(F.idr)+1),                ...
            'The number of nodes and of elements don''t match' ...
            );
         
         return
         
         % check the numbering of DOFs (if any)
         t_ = [];
         for v = 1 : numel(F.Dofs)
            t_ = union(t_, F.Dofs{v});
         end
         assert(max(t_(:)) == F.nbDofs);     % The total nb of DOFs matches max(Dofs)
         assert(all(diff(t_) == 1));         % Numbering of DOFs is consecutive
      end
   end
   
   methods(Static)
      function id = getNextFreeID(listOf_id)
         % getNextFreeID               Return the next free ID from a liste of IDs (which may eventually be empty)
         %
         % Description:
         %
         % Input:
         % - listOf_id                 List of IDs which already exist (if any)
         %
         % Output:
         % - id                        Next free ID
         %
         % Notes:                  1)  This function is for internal usage only. It is not at all optimized.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 09-Apr-2017 - first version.
         
         if isempty(listOf_id)
            id = 1;
         else
            id = 1 + max(listOf_id(:));
         end
      end
   end
end
