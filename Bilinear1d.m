classdef Bilinear1d < handle
   
   % Bilinear1d                  Bilinear form
   %
   % Description:                This class describes a bilinear form like (a, b) which can be evaulated
   %                             on quadrature points. This class provides a simpler way for writing
   %                             FE weak formulations.
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
   % 13-Apr-2017 - First version.
   
   % --------------------------->| description of the function ---|------------------------------------------->| remarks
   
   properties(GetAccess='public', SetAccess='protected')
      op=+1;                     % operator (it may be +1 or -1)
      lhs_term;                  % left-hand-side term of the bilinear form
      rhs_term;                  % right-hand-side  "  "   "    "       "
   end
   
   methods
      function B = Bilinear1d(lhs, rhs)
         % Bilinear1d                  Constructor
         %
         % Description:                This is the constructor of bilinear forms (a,b). One, and only one,
         %                             among a and b MUST be marked as test-function.
         %
         %
         % Input:
         % - lhs                       Left-hand-side  (aka a)
         % - rhs                       Right-hand-side (aka b)
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
         % 14-Apr-2017 - first version.
         
         % This constructor allows fixing a bug, but it is very questionable
         if nargin == 0 , return ; end
         
         if nargin == 1  &&  isa(lhs, 'Bilinear1d')
            % Copy constructor
            B_ = lhs;            % = original object
            B = Bilinear1d.empty();
            for k = 1 : numel(B_)
               B(k).op = B_(k).op;
               B(k).lhs_term = B_(k).lhs_term;
               B(k).rhs_term = B_(k).rhs_term;
            end
            return
         end
         
         nb = 0;
         if isa(lhs, 'ShapeFun1d')  &&  lhs.testFun , nb = nb + 1 ; end
         if isa(rhs, 'ShapeFun1d')  &&  rhs.testFun , nb = nb + 1 ; end
         assert(nb == 1, 'One (and only one) of the terms of a bilinear form MUST be a test-function');
         
         if isa(lhs, 'ShapeFun1d') &&  lhs.testFun
            t_ = lhs ; lhs = rhs ; rhs = t_;
         end
         
         B.op = +1;
         B.lhs_term = lhs;
         B.rhs_term = rhs;
      end
      
      % Algebrical operators
      function B = plus(B1, B2)
         % plus                        Plus operator
         %
         % Description:                This method returns the sum of two bilinear forms.
         %
         % Input:
         % - B1                        left-hand-side bilinear forms
         % - B2                        right  "   "      "       "
         %
         % Output:
         % - B                         Bilinear form B1 + B2
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
         % 14-Apr-2017 - first version.
         
         % B = Bilinear1d.empty();
         % for k = 1 :numel(B1)
         %    B(end+1) = Bilinear1d(B1(k));
         % %end
         B = Bilinear1d(B1);
         for k = 1 :numel(B2)
            B(end+1) = Bilinear1d(B2(k));
         end
      end
      function B = uplus(B1)
         % uplus                       Unitary plus operator
         %
         % Description:
         %
         % Input:
         % - B1                        Bilinear form
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
         % 14-Apr-2017 - first version.
         
         B = Bilinear1d(B1);
      end
      function B = uminus(B1)
         % uminus                      Unitary minus operator
         %
         % Description:
         %
         % Input:
         % - B1                        Bilinear form
         %
         % Output:
         % - B                         Opposite of the bilinear form
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
         % 14-Apr-2017 - first version.
         
         B = Bilinear1d(B1);
         for k = 1 : numel(B)
            B(k).op = - B(k).op;
         end
      end
      function B = minus(B1, B2)
         % minus                       Minus operator
         %
         % Description:                This method computes the bilinear form:  A + (-B)
         %
         % Input:
         % - B1                        left-hand-side bilinear forms
         % - B2                        right  "   "      "       "
         %
         % Output:
         % - B                         Bilinear form B1 - B2 = B1 + (-B2)
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
         % 14-Apr-2017 - first version.
         
         B = B1 + (-B2);
      end
      
      function [dofs_i, dofs_j] = getDofs(B, FEM)
         % getDofs                     Get the two tables of DOFs of a bilinear form (a,b)
         %
         % Description:
         %
         % Input:
         % - B                         Bilinear form
         % - FEM                       FEm
         %
         % Output:
         % - dofs_i                    Table of DOFs of the l.h.s. term (if any)
         % - dofs_j                      "   "   "   "   "  r.h.s. term (= test function)
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
         % 14-Apr-2017 - first version.
         if isa(B.lhs_term, 'ShapeFun1d')
            dofs_i = B.lhs_term.getDofs(FEM);
         else
            dofs_i = [];
         end
         dofs_j = B.rhs_term.getDofs(FEM);
      end
      function [A, b] = feval(B, varargin)
         % feval                       Evaluate a single bilinear form
         %
         % Description:
         %
         % Input:
         % - B                         Bilinear form
         % - varargin                  Additional arguments to pass to bilinear@ShapeFun1d
         %
         % Output:
         % - A(e,q,i,j)                Element-wise stiffness matris
         % - b(e,q,1,j)                   "     "   force vector
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
         % 14-Apr-2017 - first version.
         
         [A, b] = bilinear(B.lhs_term, B.rhs_term, varargin{:});
         assert(B.op == 1  |  B.op == -1);
         if B.op == -1 , A = -A ; b = -b ; end
      end
   end
end
