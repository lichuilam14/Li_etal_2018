
function y = vec_discrete(x,s,n)
    %% [y] = VEC_DISCRETE(x,s,n)
    % discretize vector [x] in [n] bins, independently for each index in [s]
    % if [s] is empty, the discretization is processed globally
    
    %% function
    
    % default
    func_default('s',ones(size(x)));
    
    % assert
    assertVector(x);
    assertVector(s);
    assertScalar(n);
    assertScalar(n);
    assertSize(x,s);
    
    % do
    y = nan(size(x));
    [u_s,n_s] = numbers(s);
    for i_s = 1:n_s
        ii_s = (s == u_s(i_s));
        l_s  = sum(ii_s);
        x_s  = x(ii_s);
        [~,ss_s] = sort(x_s(:));
        y_s = nan(1,l_s);
        y_s(ss_s) = ceil((1:l_s) * n / l_s);
        y(ii_s) = y_s;
    end
end
