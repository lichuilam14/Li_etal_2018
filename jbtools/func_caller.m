
function c = func_caller(n)
    %% c = FUNC_CALLER([n])
    % returns the name of the caller function
    % n : number of levels ahead (default 1)
    % c : string with the name of the caller
    
    %% function
    
    % default
    if ~nargin, n = 1; end
    n = n+2; % correction
    
    % assert
    assert(n>=0,'func_caller: error. n must be positive');
    
    % get name
    db = dbstack();
    if numel(db)<n, c = '';
    else            c = db(n).name;
    end
    
end