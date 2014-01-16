function x = apply(f,varargin)
% If f is a function that takes two inputs and returns an output of the 
% same type, this repeatedly applies said function across an arbitrary 
% number of inputs, from right to left. If the function is commutative and 
% associative (i.e. addition, multiplication, intersection, union) the
% result will not depend on the order.

if length(varargin) > get(0,'RecursionLimit')
    set(0,'RecursionLimit',length(varargin));
end
% It would be really cool if MATLAB could do tail recursion so this could
% run in constant space and still be implemented recursively.

switch length(varargin)
    case 0
        x = [];
    case 1
        x = varargin{1};        
    otherwise
        x = f(varargin{1},apply(f,varargin{2:end})); % Basically all I want in life is for MATLAB to actually be LISP.
end