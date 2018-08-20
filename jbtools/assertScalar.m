
function assertScalar(varargin)
    %% ASSERTSCALAR(x1,x2,..)
    % assert that variables [x#] are scalars
    
    %% function
    b = cellfun(@isscalar,varargin);
    b = logical(b);
    b = all(b);
    c = func_caller();
    func_default('c',func_caller(0));
    assert(b,'%s: error. one or more variables are not scalar',c);

end