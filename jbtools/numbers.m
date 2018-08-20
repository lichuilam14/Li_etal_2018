
function [u,n,h] = numbers(s)
    %% [u,n,h] = NUMBERS(s)
    % s can be a single vector or a struct
    % if a struct, numbers will be applied to all fields within s

    %% function
    % single vector
    if ~isstruct(s)
        [u,n,h] = get_unh(s,nargout);
        return;
    end
    
    % struct
    u = struct();
    n = struct();
    h = struct();
    f = fieldnames(s);
    for i = 1:length(f)
        [u.(f{i}),n.(f{i}),h.(f{i})] = get_unh(s.(f{i}),nargout);
    end
end

%% auxiliar
function [u,n,h] = get_unh(x,k)
    if ~iscell(x), x = vec_filter(x); end
    u = unique(x);
    n = length(u);
    h = nan(n,1);
    if k == 3
        if ~iscell(x),
            for i=1:n, h(i) = sum(x==u(i)); end
        else
            for i=1:n, h(i) = sum(strcmp(x,u{i})); end
        end
    end
end
