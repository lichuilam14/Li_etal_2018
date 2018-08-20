
function assertSize(varargin)
    %% ASSERTSIZE(x1,x2,..)
    % assert that all arguments have the same size
    
    %% function
    c = func_caller();
    if isempty(c), c = 'assertSize'; end % func_default('c','assertSize'); ... infinite loop
    s1 = size(varargin{1});
    n1 = length(s1);
    for i = 2:nargin
        si = size(varargin{i});
        ni = length(si);
        assert(n1==ni && all(s1==si),'%s: error. one or more variables have different size',c);
    end
end
