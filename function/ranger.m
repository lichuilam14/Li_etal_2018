
function r = ranger(x)
    %% r = RANGER(x)
    
    %% function
    x = mat2vec(x);
    r_min = nanmin(x);
    r_max = nanmax(x);
    if isempty(r_min), r_min = nan; end
    if isempty(r_max), r_max = nan; end
    r = [r_min,r_max];
end
