
function y = vec_filter(x,f)
    %% [y] = VEC_FILTER(x,f)
    % filter elements from a vector
    % f : filter function (default "~isnan")
    % x : input vector
    % y : output (sub-)vector
    % example
    %   >> y = vec_filter(x,@(a) a>0);
    
    %% function
    func_default('f',@(a)~isnan(a));
    y = x(f(x));
end
