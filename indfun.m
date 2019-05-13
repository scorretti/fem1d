function f = indfun(x, A)

% indfun                      Indicator function
% 
% Description:                The indicator function returns 1 for all x belonging to A, and 0 otherwise.
%                             In the case of this function, A = [a b] is considered, therefore:
%
%                                f = { 1     if a <= x <= b, 
%                                      0     otherwhise
%
% Input:
% - x                         x
% - A                         Interval in the form [a b]
%
% Output:                     
%  -f                         1 if x belongs to [a, b], 0 otherwise
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

a = A(1) ; b = A(2);
f = zeros(size(x));
f(a <= x  &  x <= b) = 1;

end


