classdef ShapeFun1d < handle
   
   % ShapeFun1d                  This class implements 1D shape functions
   %
   % Description:                In 1D, shape functions may be associated with the nodes, or with the
   %                             single edge of each element. The properties fn and fa are cell-array
   %                             which contain the handle to shape functions associated respectively with
   %                             the n-th node (fn{n,k} with n = 1 or 2) and with the edge (fa{1,k}).
   %
   % Properties:
   % - fn{n,k}                   Handle to the k-th shape function associated with the n-th node
   % - fa{1,k}                      "   "   "  "      "       "        "       "   the only edge
   % - d_ord                     Order of differentiation
   % - prefa                     Prefactor which accounts for the sign (it is allowed to be +1 or -1)
   % - testFun                   Flag: true if the shape function is used as test function
   % - Fem                       FEMs linked with this variable
   % - id                        ID-number of this variable with respect of each FEMs
   %
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
   % 10-Apr-2017 - Prefactor introduced.
   % 13-Apr-2017 - Modification: the class Bilinear1d is taken into account.
   % 13-Apr-2017 - Modification: a link between Fem1d and ShapeFun1d is created.
   % 16-Apr-2017 - Modification: the jacobian is retrieved directly from the FEM.
   
   % --------------------------->| description of the function ---|------------------------------------------->| remarks
   
   properties(GetAccess='public', SetAccess='protected')
      fn = {};                   % = handle to functions associates to nodes
      fa = {};                   % = handle to functions associates to nodes
      d_ord = 0;                 % = order of derivation
      prefa = 1;                 % = prefactor
      testFun = false;           % = flag: true if the shape function is used as test function
      Fem = Fem1d.empty();       % = Fem which are associated with this shape function
      id = [];                   % = ID-numbers of this shape function with respect of each Fem
   end
   
   methods
      function S = ShapeFun1d(varargin)
         % ShapeFun1d                  Constructors
         %
         % Description:
         %
         % Input: (copy constructor)
         % - S_                        Shape functions to be copied
         %
         % Output:
         % - S                         Deep copy of the shape functions
         %
         %
         % Input: (Lagrange functions of any order)
         % - 'P'
         % - order                     Order of Lagrange functions
         %
         % Input: (Lagrange functions of any order)
         % - 'P1', 'P2' ...            Same as the previous synapsys
         %
         % Input: (custom shape functions)
         % - 'custom'
         % - sn{}                      Shape functions associated with nodes
         % - sa{}                        "       "        "        "   edges
         %
         % Output:
         % - S                         Shape functions
         %
         %
         % Input: (for internal usage only)
         % - '*', or 'dummy'
         %
         % Output:
         % - S                         Dummy test functions of type P1 (for internal usage only)
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
         % 10-Apr-2017 - improvement:  the prefactor is implemented (in the copy constructor).
         % 13-Apr-2017 - Modification: a link between Fem1d and ShapeFun1d is created.
         
         if nargin == 1  &&  isa(varargin{1}, 'ShapeFun1d')
            % Copy constructor
            S_ = varargin{1};          % = original object
            for k = 1 : numel(S_)
               S(k).fn = S_(k).fn;
               S(k).fa = S_(k).fa;
               S(k).d_ord = S_(k).d_ord;
               S(k).prefa = S_(k).prefa;
               S(k).testFun = S_(k).testFun;
               S(k).Fem = S_(k).Fem;
               S(k).id = S_(k).id;
            end
            return
         end
         
         % This is a shortcut which adds no new capabilities to the class
         if nargin == 1  &&  ischar(varargin{1})
            t_ = varargin{1};
            type = upper(t_(1));
            order = str2double(t_(2:end));
            S = ShapeFun1d(type, order);
            return
         end
         
         % Main switch
         type = upper(varargin{1});
         switch type
            case 'P'                   % Lagrange shape functions
               order = varargin{2};    % = order of the shape functions
               
               % Handle the special case P0
               if order == 0
                  S.fn = {};
                  S.fa = {@(u) 1};
                  return;
               end
               
               f_ = {};
               u = linspace(0, 1, 1+order);     % = coordinates of the nodes in the reference space
               
               % Compute Lagrange shape functions for each node
               for k = 1 : order+1
                  fk = zeros(size(u)) ; fk(k) = 1;
                  coeffs = polyfit(u, fk, order);
                  f_{k} = @(u) polyval(coeffs, u);
               end
               
               % Assign the shape functions to the properties fn and fa
               S.fn{1,1} = f_{1} ; S.fn{2,1} = f_{end};
               S.fa = f_(2:end-1);
               
            case {'*', 'dummy'}        % dummy test function (for internal usage only)
               S = ShapeFun1d('P0');
               S.testFun = true;
               
            case 'CUSTOM'              % User-provided custom shape functions
               S.fn = varargin{2};
               if nargin >= 3 , S.fa = varargin{3} ; end
               
            otherwise
               error('Unknown type');
         end
      end
      
      % Association between FEM and shape functions
      function associate(S, FEM, id)
         % associate                   Create an association between the shape function and a FEM
         %
         % Description:                This method creates an link between this shape function and a FEM.
         %                             This method, for internal usage only, makes the shape funcion
         %                             "aware" of the fact that is is being used as unknown for a given
         %                             FEM. Note that the table of DOFs is stored with the FEM. This means
         %                             that a single shape function can be associated with different FEMs,
         %                             with a different table of DOFs for each FEM.
         %
         %                             This method is intended for internal usage only.
         %
         %
         % Input:
         % - S                         Shape function
         % - FEM                       FEM
         % - id                        ID-number of this shape function with respect of FEM
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
         
         % --------------------------->| description -------------------|------------------------------------------->| remarks
         
         t_ = S.getIdOfFem(FEM);
         if t_ ~= 0
            warning(sprintf([                                                       ...
               'The shape function is already associated with the FEM\n'            ...
               '         (this operation will override previous associations)']));
         else
            S.Fem(end+1) = FEM;
            t_ = numel(S.Fem);
         end
         S.id(t_) = id;
      end
      function id = getIdOfFem(S, FEM)
         % getIdOfFem                  ID-number of the FEM (for internal usage only)
         %
         % Description:                This method returns the ID-number of the provided FEM.
         %
         %                             If this shape function is being used as unknown for a given FEM
         %                             (that is, if it is associated with the FEM), then a copy of the
         %                             handler to that FEM is stored in the property Fem. In this case,
         %                             this method will return a positive index. Otherwise, this method
         %                             will return 0.
         %
         %
         % Input:
         % - S                         Shape function
         % - FEM                       FEM
         %
         % Output:
         % - id                        ID-number of the FEM with respect of the shape function
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
         
         [~, id] = ismember(FEM, S.Fem);
      end
      function id = getId(S, FEM)
         % getId                       ID of the shape function with respect of a FEM (for internal usage)
         %
         % Description:                This method returns the ID-number of this shape function with
         %                             respect of the given FEM, or 0 if this shape function is not
         %                             associated with the FEM.
         %
         % Input:
         % - S                         Shape function
         % - FEM                       FEM
         %
         % Output:
         % - id                        ID-number of the shape function with respect of the FEM
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
         
         t_ = S.getIdOfFem(FEM);
         if t_ == 0
            id = 0;
         else
            id = S.id(t_);
         end
      end
      
      % Operators
      function dS = d(S)
         % d                           Derivative operator
         %
         % Description:                This method exercutes a deep copy of the shape function (*) and
         %                             increase the value of the property d_ord, meaning that the value
         %                             which the method feval() must compute is the derivative of the
         %                             shape functions, and not the shape functions theirselves.
         %
         %                         (*) It is mandatory to execute a deep copy because some of the
         %                             properties are modified.
         %
         % Input:
         % - S                         Shape functions
         %
         % Output:
         % - dS                        Derivative of the shape functions
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
         
         dS = ShapeFun1d(S);
         dS.d_ord = dS.d_ord + 1;
      end
      function aS = mtimes(a, S)
         % mtimes                      Multiply a shape function by a prefactor
         %
         % Description:                This method has in fact two different synosys:
         %
         %                             The first one, returns the shape function resulting of the product of
         %                             a prefactor by a shape function.
         %                             In practice, this method makes a deep copy of the shape function
         %                             and modifies the property prefa of the copy.
         %
         %                             The second one, takes a constant, a function handler or a shape
         %                             function AND another shape function (used as test function) and
         %                             returns a bilinear form.
         %
         %                             Please note that Fem1d is a simple "toy" for teaching FEM, and not
         %                             a general purpose FEM environment. Therefore its capabilites are
         %                             minimalistic. In particuliar, it is possible to multiply only
         %                             one time a shape function by a prefactor.
         %
         % Input:
         % - a                         Prefactor
         % - S                         Shape functiàon
         %
         % Output:
         % - NaS                       Product of prefactor * shape function
         %
         %
         % Input:
         % - S1                        l.h.s. term of the bilinear form
         % - S2                        r.h.s. term of the bilinear form
         %
         % Output:
         % - B                         Bilinear form
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
         % 13-Apr-2017 - Modification: the class Bilinear1d is taken into account.
         
         if (isa(a, 'ShapeFun1d')  &&  a.testFun)  ||  (isa(S, 'ShapeFun1d')  &&  S.testFun)
            aS = Bilinear1d(a, S);
         else
            aS = ShapeFun1d(S);
            aS.prefa = a;
         end
      end
      function Sp = ctranspose(S)
         % ctranspose                  Mark the shape function as test function
         %
         % Description:                This method makes a deep copy of the shape function, and marks it
         %                             as shape function (property testFun = true). When creating a new
         %                             bilinear form, at least one argument MUST be a shape function
         %                             marked as test function.
         %
         % Input:
         % - S                         Shape function
         %
         % Output:
         % - Sp                        Shape function (marked as test function)
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
         
         Sp = ShapeFun1d(S);
         Sp.testFun = true;
      end
      
      % These operators are only useful to define constraints
      function B = plus(S1, S2)
         % plus                        Plus operator
         %
         % Description:                This method takes two terms (one of which must be a shape function)
         %                             and returns the dummy bilinear form:
         %
         %                                (S1, *) + (S2, *)
         %
         %                             The only purpose of this and other operators is to allow writing
         %                             the constraints in a more user-friendly formalism.
         %
         % Input:
         % - S1                        First term
         % - S2                        Second term
         %
         % Output:
         % - B                         Dummy bilinear form (S1, *) + (S2, *)
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:                   impose@Fem1d
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         B = S1*ShapeFun1d('*') + S2*ShapeFun1d('*');
      end
      function B = minus(S1, S2)
         % minus                       Minus operator
         %
         % Description:                This method takes two terms (one of which must be a shape function)
         %                             and returns the dummy bilinear form:
         %
         %                                (S1, *) - (S2, *)
         %
         %                             The only purpose of this and other operators is to allow writing
         %                             the constraints in a more user-friendly formalism.
         %
         % Input:
         % - S1                        First term
         % - S2                        Second term
         %
         % Output:
         % - B                         Dummy bilinear form (S1, *) - (S2, *)
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:                   impose@Fem1d
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         B = S1*ShapeFun1d('*') - S2*ShapeFun1d('*');
      end
      function B = uplus(S1)
         % uplus                       Unitary plus operator
         %
         % Description:                This method takes a shape function and returns the dummy bilinear form:
         %
         %                                (S1, *)
         %
         %                             The only purpose of this and other operators is to allow writing
         %                             the constraints in a more user-friendly formalism.
         %
         % Input:
         % - S1                        First term
         %
         % Output:
         % - B                         Dummy bilinear form (S1, *)
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:                   impose@Fem1d
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         B = S1*ShapeFun1d('*');
      end
      function B = uminus(S1)
         % uminus                      Unitary minus operator
         %
         % Description:                This method takes a shape function and returns the dummy bilinear form:
         %
         %                                (-S1, *)
         %
         %                             The only purpose of this and other operators is to allow writing
         %                             the constraints in a more user-friendly formalism.
         %
         % Input:
         % - S1                        First term
         %
         % Output:
         % - B                         Dummy bilinear form (-S1, *)
         %
         % Notes:
         %
         % Exemple:
         %
         % See also:                   impose@Fem1d
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         B = -(S1*ShapeFun1d('*'));
      end
      function B = eq(S1, S2)
         % eq                          Equality operator
         %
         % Description:                This method takes two terms (one of which must be a shape function)
         %                             and returns the dummy bilinear form:
         %
         %                                (S1, *) - (S2, *)
         %
         %                             The only purpose of this and other operators is to allow writing
         %                             the constraints in a more user-friendly formalism.
         %
         % Input:
         % - S1                        First term
         % - S2                        Second term
         %
         % Output:
         % - B                         Dummy bilinear form (S1, *) - (S2, *)
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         if isa(S2, 'ShapeFun1d')
            B = S1 - S2;
         else
            % Special case (to generalize)
            B = S1*ShapeFun1d('*') - Bilinear1d(S2,ShapeFun1d('*'));
         end
      end
      
      % Methods which deals with DOFs
      function [Nn, Na] = getMultiplicity(S)
         % getMultiplicity             Return the number of DOFs on each node and edge
         %
         % Description:
         %
         % Input:
         % - S                         Shape functions
         %
         % Output:
         % - Nn                        Number of DOFs associated with each node
         % - Na                           "   "   "       "       "   each edge
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
         
         Nn = size(S.fn, 2);        % = nb of shape functions for each node
         Na = numel(S.fa);          % = nb of shape functions for each edge
      end
      function [dofs, nb] = buildDofs(S, FEM, ide)
         % buildDofs                   Generate the DOFs for the given shape functions
         %
         % Description:                This method generates the DOFs for a subset of elements of a 1D mesh.
         %                             For each element, DOFs are numbered in the following order:
         %                              - DOFs on node 1,
         %                              - DOFs on node 2,
         %                              - DOFs on the edge.
         %
         % Example:                    For instance, consider the following mesh composed of 3 elements:
         %
         %                                         0 ------- 1 ------- 2 ------- 3
         %
         %                             For P1 functions (1 DOF on each node, no DOFs on edges) on has:
         %
         %                                         0 ------- 1 ------- 2 ------- 3
         %                                        (1)       (2)       (3)       (4)
         %
         %                             >> FEM = Fem1d(0:3) ; S = ShapeFun1d('P1');
         %                             >> dofs = S.buildDofs(FEM, 1:3)
         %
         %                             dofs =
         %
         %                                  1     2      <--- DOFs on the first element
         %                                  2     3
         %                                  3     4
         %
         %                             On the same mesh, for P4 functions (1 DOF on each node, 3 DOFs on
         %                             each edge) on has:
         %
         %                             >> FEM = Fem1d(0:3) ; S = ShapeFun1d('P1');
         %                             >> dofs = S.buildDofs(FEM, 1:3)
         %
         %                             dofs =
         %
         %                                  1     2     5     6     7      <--- DOFs on the first element
         %                                  2     3     8     9    10
         %                                  3     4    11    12    13
         %
         %                                             |<----------->|     = DOFs on the edge
         %                               |<->|  |<->|                      = DOFs on the nodes
         %
         %
         %
         % Input:
         % - S                         Shape function
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       ID-number of the elements for which DOFs have to be generated
         %
         % Output:
         % - dofs(e,p)                 DOFs associated with each element
         % - nb                        Number of DOFs
         %
         %
         % Notes:                  1)  In fact, the argument FEM is not used.
         %
         %                         2)  Only the case of a single DOF on each node has been tested.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 12-Apr-2017 - first version.
         
         
         assert(isa(FEM, 'Fem1d')  &&  (max(ide) <= FEM.getNbOfElements())  &&  (min(ide) >= 0));
         
         [Nn, Na] = S.getMultiplicity();
         Ne = numel(ide);
         dofs = zeros(Ne,2*Nn+Na);
         
         % Generate DOFs on nodes
         for k = 1 : 2*Nn
            dofs(:,k) = k + Nn*(0 : Ne-1)';
         end
         
         % Generate DOFs on edges
         t_ = 2*Nn+Nn*(Ne-1) ; assert(t_ == max(dofs(:)));
         for k = 1 : Na
            dofs(:,2*Nn+k) = t_ + k + Na*(0 : Ne-1)';
         end
         
         % Compute the number of  DOFs
         nb = 2*Nn+Nn*(Ne-1) + Na + Na*(Ne-1) ; assert(nb == max(dofs(:)));
      end
      function dofs = getDofs(S, FEM)
         % getDofs                     Get the table of DOFs of the shape function
         %
         % Description:                This method returns the table of DOFs of the shape function with
         %                             respect of the given FEM. If the shape funcion is not associated
         %                             with that FEM, this method will return [].
         %
         %                             The table of DOFs (dofs(e,p)) is provided for all of the elements.
         %                             If the shape function is not declared on some elements, the
         %                             corresponding entries in the table of DOFs will be 0.
         %
         % Input:
         % - S                         Shape function
         % - FEM                       FEM
         %
         % Output:
         % - dofs(e,p)                 Table of DOFs
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
         
         t_ = S.getIdOfFem(FEM);       % = ID-number of the FEM
         
         if t_ > 0
            % The shape function is associated with the FEM
            id = S.id(t_);
            dofs = S.Fem(t_).getDofs(id);
            
         else
            % The shape function is not associated with the FEM
            dofs = [];
         end
      end
      
      % Methods which perform some computations
      function [vn, va] = feval(S, FEM, ide, u, solution)
         % feval                       Evaluate the shape functions in the quadrature points
         %
         % Description:                This method evaluates the shape functions, or their derivatives
         %                             (according to the value of the property d_ord) on some quadrature
         %                             points, the coordinates of which are provided with respect of the
         %                             reference element, and therefore are the same for each element of
         %                             the real space.
         %
         %                             The method accepts as input a Finite Element discretization (FEM),
         %                             a list of ID-number of elements (ide) and the coordinates u(q)
         %                             in the reference space (q = index of the quadrature points).
         %
         %                             Optionally, a jacobian may be provided. By default, the standard
         %                             linear transformation that maps [0 1] --> [Xn(1) Xn(2)] will be used.
         %
         %                             The method returns two matrix:
         %                              - vn(e,q,n,k) = value of the shape functions computed in the q-th
         %                                quadrature point of the e-th element, for the k-th shape
         %                                function associated with the n-th node.
         %                              - va(e,q,1,k) = value of the shape functions computed in the q-th
         %                                quadrature point of the e-th element, for the k-th shape
         %                                function associated with the only edge of the element.
         %
         %                             If a single output argument is provided, the method returns a
         %                             matrix vna(e,q,p,k). In this case, the order of the DOFs is the
         %                             same as explained in the method buildDofs(), that is:
         %                              - DOFs on node 1,
         %                              - DOFs on node 2,
         %                              - DOFs on the edge.
         %
         %  Example:                   For instance, consider the following mesh:
         %
         %                                         0 ------- 1 ------- 2 ------- 3
         %
         %                             >> FEM = Fem1d(0:3) ; S = ShapeFun1d('P3');
         %
         %                             For a given value of u, the shape functions will take the same value
         %                             on each element. Consider the case u = 0, which corresponds to the
         %                             following situation:
         %
         %                                S  (0) = 1     ;     S  (0) = S  (0) = S  (0) = 0
         %                                 n1                   n2       a1       a2
         %
         %                             One has:
         %
         %                             >> vna = S.feval(FEM, ide, 0) ; reshape(vna, [3 4])
         %
         %                             ans =
         %
         %                                 1.0000    0.0000    0.0000   -0.0000
         %                                 1.0000    0.0000    0.0000   -0.0000
         %                                 1.0000    0.0000    0.0000   -0.0000
         %
         %                                 |<-->|                                = S  (u)
         %                                                                          n1
         %
         %                             For u = 1/3 on expects that:
         %
         %                                S  (0) = 1     ;     S  (0) = S  (0) = S  (0) = 0
         %                                 e1                   n1       a2       a2
         %
         %                             Indeed:
         %
         %                             >> vna = S.feval(FEM, ide, 0) ; reshape(vna, [3 4])
         %
         %                             ans =
         %
         %                                 0   -0.0000    1.0000   -0.0000
         %                                 0   -0.0000    1.0000   -0.0000
         %                                 0   -0.0000    1.0000   -0.0000
         %
         %                                                |<-->|              	= S  (u)
         %                                                                          a1
         %
         % Interpolation:              If a vector containing the solution for the linear system is
         %                             provided, this method dispatch the execution to feval2() which 
         %                             computes the solution in the given points (u) on each element (ide)
         %
         %
         % Input:
         % - S                         Shape functions
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       Id-number of the nodes where the shape functions have to be evaulated
         % - u                         Coordinates of the quadrature points (can be a vector)
         %
         % Output:
         % - vn(e,q,n,k)               Value of the shape functions associated with the nodes
         % - va (e,q,1,k)                "   "   "    "       "          "      "    "  edges
         %
         % Output:
         % -vna(e,q,p)                 Value of the shape functions associated with nodes and edges
         %                             "packed" together
         %
         %
         % Input:
         % - S                         Shape functions
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       Id-number of the nodes where the shape functions have to be evaulated
         % - u                         Coordinates of the quadrature points (can be a vector)
         % - solution                  Solution of the linear system A.x = b
         %
         % Output:
         % -val(e,q)                   Value of the expanded quantity on each computation point
         % - x(e,q)                    Coordinates of the computation points
         % 
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
         
         % If the values of DOFS are provided, dispatch the execution to feval2()
         if nargin == 5  &&  ~ isempty(solution)
            [vn, va] = S.feval2(FEM, ide, u, solution);
            return
         end
         
         J = FEM.getJacobian(ide);              % = jacobian which maps: [0,1] --> real elements
         
         % Pre-allocate the variables vn, ve
         Nq = numel(u);                         % = nb of computation (quadrature) points for each element
         [Nn, Na] = S.getMultiplicity();        % = nb of shape functions for each node (Nn), edge (Na)
         Ne = numel(ide);                       % = nb of elements
         
         vn = zeros(Ne, Nq, 2, Nn);             % = value of vn in the e-element, q-th point, k-th function
         va = zeros(Ne, Nq, 1, Na);             % =   "   "  ve "   "      "      "     "     "    "
         [~ , Xn] = FEM.getElements('id', ide); % Xn(e,n) = coordinates of the nodes in the real space
         
         for e = 1 : Ne
            % Compute the shape functions associated with the nodes
            for n = 1 : 2
               for k = 1 : Nn
                  vn(e,:,n,k) = S.feval_(S.prefa, S.fn{n,k}, S.d_ord, u, Xn(e,:), J);
               end
            end
            
            % Compute the shape functions associated with the only edge (that's why se set a = 1)
            a = 1;
            for k = 1 : Na
               va(e,:,a,k) = S.feval_(S.prefa, S.fa{a,k}, S.d_ord, u, Xn(e,:), J);
            end
         end
         
         % If a single output argument is required, cat together vn and va
         if nargout == 1
            vn = vn(:,:,:) ; va = va(:,:,:);
            
            if isempty(vn)
               vn = va;
            elseif isempty(va)
               ;
            else
               vn = cat(3, vn, va);
            end
         end
      end
      function [A, b] = bilinear(S1, S2, FEM, ide, u)
         % bilinear                    Evaluate a bilinear form (S1, S2) in the quadrature points
         %
         % Description:                This method evaluates a bilinear form (S1, S2).
         %
         %                             The first term S1 may be a constant value, the handle to a function
         %                             f(x), or a shape function.
         %                             The second term Sé must be a shape function, which plays the role
         %                             of test function.
         %
         %                             If the first term is a shape function S1(x), then this method will
         %                             returns b = [] and a matrix A(e,q,i,j):
         %
         %                                A(e,q,i,j) = (S1 (x) , S2 (x))
         %                                                i        j
         %
         %                             where: e = index of elements, q = index of quadrature points.
         %
         %                             If the first term is a source term (i.e. a constant value K or the
         %                             handle to a function f(x)), this method will returns A = [] and
         %                             a matrix b(e,q,1,j):
         %
         %                                b(e,q,1,j) = (K , S2 (x))
         %                                                    j
         %                             or
         %
         %                                b(e,q,1,j) = (f(x) , S2 (x))
         %                                                       j
         %
         %                             where: e = index of elements, q = index of quadrature points.
         %
         %                             In order to finalize the resolution of a Finite Elelent problem,
         %                             the returned matrix must be "integrated" so as to remove the q
         %                             dimension:
         %
         %                                A(e,q,i,j)  -->  A(e,1,i,j)
         %                                b(e,q,1,j)  -->  b(e,1,1,j)
         %
         %                             and then assembled respectively in the stiffness matrix (as for A)
         %                             ond in the force vector (as for b).
         %
         %                             If in place of the coordinates u the name of a quadrature formula
         %                             is provided, this method returns directly the integrated matrix.
         %
         %
         % Input:
         % - S1                        Left-hand-side term  (constant, handle_function or ShapeFun1d)
         % - S2                        Right-hand-side term (ShapeFun1d)
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       Id-number of the nodes where the bilinear form has to be evaulated
         % - u                         Coordinates of the quadrature points (can be a vector)
         %
         % Output:
         % - A(e,q,i,j)                Element-wise stiffness matrix
         % - b(e,q,1,j)                Element-wise force vector
         %
         %
         %
         % Input:
         % - S1                        Left-hand-side term  (constant, handle_function or ShapeFun1d)
         % - S2                        Right-hand-side term (ShapeFun1d)
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       Id-number of the nodes where the bilinear form has to be evaulated
         % - quad                      Name of a quadrature formula ('QUAD3', 'QUAD4' ... 'QUAD6', 'QUAD10'), or
         %                             quadrature formula {u, w}, where: u = coordinates of Gauss points,
         %                             w = weights
         % - J                         Jacobian (default = linear transformation)
         %
         % Output:
         % - A(e,1,i,j)                Element-wise integrated stiffness matrix
         % - b(e,1,1,j)                Element-wise integrated force vector
         %
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
         % 12-Apr-2017 - first version.
         % 17-Apr-2017 - Improvement:  more quadrature formulas added (up to QUAD6).
         
         J = FEM.getJacobian(ide);              % = jacobian which maps: [0,1] --> real elements
         
         % If a quadrature formula is provided, compute directly the integrated matrix
         if ~isnumeric(u)
            if ischar(u)
               % The name of a quadrature formula is provided
               switch lower(u)
                  case 'quad3'
                     u = [ 0.112701665379258   0.500000000000000   0.887298334620742 ];
                     w = [ 0.277777777777778   0.444444444444444   0.277777777777778 ];
                  case 'quad4'
                     u = [ 0.069432000000000   0.330009500000000   0.669990500000000   0.930568000000000 ];
                     w = [ 0.173927500000000   0.326072500000000   0.326072500000000   0.173927500000000 ];
                  case 'quad5'
                     u = [ 0.953089922969332   0.769234655052841   0.500000000000000   0.230765344947158  ...
                        0.046910077030668 ];
                     w = [ 0.118463442528095   0.239314335249683   0.284444444444444   0.239314335249683  ...
                        0.118463442528095 ];
                  case 'quad6'
                     u = [ 0.966234757101576   0.830604693233132   0.619309593041598   0.380690406958402  ...
                        0.169395306766868   0.033765242898424 ];
                     w = [ 0.085662246189585   0.180380786524069   0.233956967286345   0.233956967286345  ...
                        0.180380786524069 0.085662246189585 ];
                  case 'quad10'
                     u = [ 0.986953264258586   0.932531683344492   0.839704784149512   0.716697697064624  ...
                        0.574437169490816  0.425562830509184   0.283302302935376   0.160295215850488      ...
                        0.067468316655508   0.013046735741414 ];
                     w = [ 0.033335672154344   0.074725674575290   0.109543181257991   0.134633359654998  ...
                        0.147762112357376  0.147762112357376   0.134633359654998   0.109543181257991      ...
                        0.074725674575290   0.033335672154344 ];
                  otherwise
                     error('Unimplemented feature');
               end
               
            elseif iscell(u)
               % A quadrature formula is provided
               t_ = u ; u = t_{1} ; w = t_{2};
            end
            
            % Compute A, b in the quadrature points
            [A, b] = bilinear(S1, S2, FEM, ide, u);
            
            % Evaluate the integrals
            dxdu = J.dxdu(u, {FEM, ide});             % = value of dx/du
            A = ShapeFun1d.quad(A, w, dxdu);
            b = ShapeFun1d.quad(b, w, dxdu);
            return
         end
         
         % Evaluate the right-hand-side S2 (= test function)
         val_S2 = S2.feval(FEM, ide, u);
         [Ne2, Nq2, Ns2] = size(val_S2);
         
         % Evaluate the left-hand-side (= test function or source term)
         if isa(S1, 'ShapeFun1d')
            % lhs is a test function
            val_S1 = S1.feval(FEM, ide, u);
            
            % Check that S1 and S2 have consistent sizes
            [Ne1, Nq1, Ns1] = size(val_S1);
            assert(Ne1 == Ne2  &&  Nq1 == Nq2);
            
            A = zeros(Ne1, Nq1, Ns1, Ns2);
            for i1 = 1 : Ns1
               for i2 = 1 : Ns2
                  A(:,:,i1,i2) = val_S1(:,:,i1) .* val_S2(:,:,i2);
               end
            end
            % b = zeros(Ne1, Nq1, 1, Ns2);
            b = [];
            
         else
            % lhs is a source term
            A = [];
            
            if isnumeric(S1)
               val_S1 = S1*ones(Ne2, Nq2, 1);
               
            elseif isa(S1, 'function_handle')
               % [~ , Xn] = FEM.getElements('id', ide); % Xn(e,n) = coordinates of the nodes in the real space
               % x = J.x(u, Xn);
               x = J.x(u, {FEM, ide});
               val_S1 = feval(S1, x);
               
            else
               error('Unimplemented feature');
            end
            
            % Check that S1 and S2 have consistent sizes
            [Ne1, Nq1, ~] = size(val_S1);
            assert(Ne1 == Ne2  &&  Nq1 == Nq2);
            
            b = zeros(Ne1, Nq1, 1, Ns2);
            for i2 = 1 : Ns2
               b(:,:,1,i2) = val_S1(:,:,1) .* val_S2(:,:,i2);
            end
         end
         
         % if nargout == 1 , A = cat(3, A, b) ; end
      end
      function [val, x] = feval2(S, FEM, ide, u, solution)
         % feval2                      Compute the interpolated quantity in the computation points
         %
         % Description:                This method is for  internal usage only. Use feval() instead.
         %
         % Input:
         % - S                         Shape functions
         % - FEM                       Finite Element model (i.e. mesh)
         % - ide                       Id-number of the nodes where the shape functions have to be evaulated
         % - u                         Coordinates of the quadrature points (can be a vector)
         % - solution                  Solution of the linear system A.x = b
         %
         % Output:
         % -val(e,q)                   Value of the expanded quantity on each computation point
         % - x(e,q)                    Coordinates of the computation points
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
         % 18-Apr-2017 - first version.
         
         % Compute the shape functions
         vna = S.feval(FEM, ide, u);      % vna(e,q,p) = value of the shape functions for each element, point
         [Ne, Nq, Np] = size(vna);        % Ne = nb of elements, Nq = nb of points, Np = nb of DOFs
         
         % Retrieve the table of DOFs, and generate the table val_u(e,p)
         dofs = S.getDofs(FEM);           % dofs(e,p) = ID of the DOFs
         
         % Keep only the DOFs on the elements (ide)
         dofs = dofs(ide,:);                                                                                % ###
         
         assert(size(dofs,1) == Ne  &   size(dofs,2) == Np);
         val_u = zeros(Ne, Np);           % = values of the DOFs (p) on elements (ide(e))
         
         % Handle the elements which take no DOF
         t_ = find(dofs > 0);
         val_u(t_) = solution(dofs(t_));
         
         % Compute the value of the expanded quantity
         val = zeros(Ne, Nq);             % val(e,q) = value in the computation points q of the element e
         for p = 1 : Np
            for q = 1 : Nq
               val(:,q) = val(:,q) + vna(:,q,p).*val_u(:,p);
            end
         end
         
         % If required, compute the coordinates of the points x(u)
         if nargin >= 2
            J = FEM.getJacobian(ide);   	% = jacobian which maps: [0,1] --> real elements
            x = J.x(u, {FEM, ide});       % = coordinates of the points
         end
      end
      
      % Render the shape functions on the reference element [0,1]
      function render(S)
         % render                      Draw the shape functions in a new figure (for didactical purposes)
         %
         % Description:                This method opens a new figure and plot all the shape functions.
         %                             SN = shape functions associated with nodes,
         %                             SE =   "       "          "      "   edges.
         %
         % Input:
         % - S                         Shape functions
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
         % 10-Apr-2017 - first version.
         
         figure
         u = linspace(0, 1, 1024);     % = coordinates in the reference space
         
         % Draw the shape functions associated with the nodes
         for k = 1 : numel(S.fn)
            plot(u, feval(S.fn{k},u), 'LineWidth', 2, 'DisplayName', sprintf('SN_%d(u)', k));
            hold on
         end
         
         % Draw the shape functions associated with the nodes
         for k = 1 : numel(S.fa)
            plot(u, feval(S.fa{k},u), 'LineWidth', 2, 'DisplayName', sprintf('SE_%d(u)', k));
         end
         box on ; grid on ; legend show
      end
      function hndl = plot(S, FEM, ide, u, solution, varargin)
         % plot                        Plot the expanded quantity over a subset of elements
         %
         % Description:                This method is intended to display the solution over a subset of
         %                             elements. The required arguments are the FEM, ID-numbers of elements
         %                             (ide) and a number of coordinates (u) with respect of the reference
         %                             element, and of course the solution.
         %
         % Input:
         % - S                         Shape function
         % - FEM                       FEM
         % - ide                       ID-numbers of elements
         % - u                         Coordinates with respect of the reference element [0, 1]
         % - solution                  Solution of the linear system A.x = b
         % - varargin                  Additional arguments to be passed to the function plot()
         %
         % Output:
         % - hndl
         %
         % Notes:                  1)  It is assumed that the elements are contiguous.
         %
         % Exemple:
         %
         % See also:
         %
         % References:
         %
         % Validation:
         %
         % 18-Apr-2017 - first version.
         
         % So as to avoid problem with discontinuous variables, enforce that 0 < u < 1
         u(u==0) = 1E-4 ; u(u==1) = 1 - 1E-4;
         
         % Compute the value of the expanded quantity val(e,q) and the coordinates of the points x(e,q)
         [val, x] = S.feval(FEM, ide, u, solution);
         
         % Sort the points with respect of x (that is why it is needed that elements are contiguous) and plot
         [x, i_] = sort(x(:)) ; val = val(i_) ; val = val(:);
         hndl = plot(x, real(val), varargin{:});
      end
   end
   
   methods(Static)
      function df = derive(f)
         % derive                      Compute a numerical approximation of the derivative of f(u)
         %
         % Description:                This method takes as input argument the handle to a function f(u)
         %                             and returns the handle to its derivative du(u), which is approximated
         %                             by a finite difference:
         %
         %                                                  f(u+h) - f(u-h)
         %                                f'(u)  ~  df(u) = ---------------
         %                                                        2h
         %
         %
         % Input:
         % - f                         Handle of the function f(u) to derive
         %
         % Output:
         % - df
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
         df = @(u) (f(u+h) - f(u-h))/(2*h);
      end
      function val = feval_(prefa, fun, k, u, Xn, J)
         % feval_                      Function for internal usage only (used by feval())
         %
         % Description:                This function computes the k-th derivative of the expression:
         %
         %                                prefa.fun(x)
         %
         %                             in the
         %                             quadrature points, the coordinates of are provided by u:
         %
         %                                 k
         %                                d
         %                                --- fun(x)
         %                                  k
         %                                dx
         %
         %
         % Input:
         % - prefa                     Prefactor
         % - fun                       Function to evaluate
         % - d_ord                     Order of the derivative
         % - u                         Coordinates of the quadrature points in the reference space
         % - Xn                        Coordinates of the nodes of the element
         % - J                         Jacobian to be used
         %
         % Output:
         % - val
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
         % 10-Apr-2017 - Modification: prefactor implemented.
         % 12-Apr-2017 - Modification: the derivative of the prefactor is NOT computed
         
         % Evaluate the prefactor and its derivatives (if required)
         if isnumeric(prefa)                       % The prefactor is a constant (the easiest case)
            val_prefa = prefa;                     % = value of the prefactor
            val_dprefa = 0;                        % = derivative of the prefactor
         elseif isa(prefa, 'function_handle')      % The prefactor is a function a(x) (the not-that-easy case)
            x = J.x(u, Xn);
            val_prefa = feval(prefa, x);
            
            % % If k > 0 the first derivative must be evaluated
            % if k > 0
            %    dprefa_ = ShapeFun1d.derive(prefa);
            %    val_dprefa = feval(dprefa_, x);
            % end
         else
            error('Unimplemented feature');
         end
         
         switch k
            case 0                     % Compute prefa.fun(x)
               val = val_prefa * feval(fun, u);
               
            case 1                     % Compute the first derivative of prefa.fun(x)
               dfun_ = ShapeFun1d.derive(fun);
               dudx = J.dudx(u, Xn);
               % val = val_prefa .* (feval(dfun_, u) .* dudx) + val_dprefa .* feval(fun, u);
               val = val_prefa .* (feval(dfun_, u) .* dudx);
               
            otherwise
               error('Unimplemented feature');
         end
      end
      function iA = quad(A, w, dxdu)
         % quad                        Numerical quadrature routine
         %
         % Description:                This method performs numerical integration of A(e,q,i,j) by using
         %                             the weights w(q):
         %                                              ___
         %                                              \
         %                                iA(e,1,i,j) = /__  w(q).A(e,q,i,j).|dxdu(e,q)|
         %                                               q
         %
         % Input:
         % - A(e,q,i,j)                Matrix to integrate
         % - w(q)                      Weights of the quadrature formula
         % - dxdu(e,q)                 Jacobian evaluated in x(e,q)
         %
         % Output:
         % - iA(e,1,i,j)               Integrated matrix
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
         % 18-Apr-2017 - first version.
         
         if isempty(A)
            iA = A;
            return
         end
         
         % Check that all the matrix have consistent size
         [Ne, Nq, Ni, Nj] = size(A);
         assert(numel(w) == Nq);
         assert(size(dxdu,1) == Ne  &&  size(dxdu,2) == Nq);
         
         % Multiply A(e,q,i,j) by the weight of the quadrature formula w(q) and accumulate into iA(e,1,i,j)
         iA = zeros(Ne, 1, Ni, Nj);
         for q = 1 : Nq
            % Multiply iA(e,1,i,j) by the absolute value of the jacobian determinant dxdu(e,q)
            for i = 1 : size(iA,3)
               for j = 1 : size(iA,4)
                  iA(:,1,i,j) = iA(:,1,i,j) + w(q) * A(:,q,i,j) .* abs(dxdu(:,q));
               end
            end
         end
      end
   end
end


