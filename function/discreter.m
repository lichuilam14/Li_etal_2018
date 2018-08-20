function out = discreter(input,varargin)
% e.g. splitter(input,10 15, 30)
bins   = [-Inf cell2mat(varargin) +Inf];
out    = ones(size(input));
level  = 1:length(bins)-1;
for b = 1:length(bins)-1
   ind = find(input >=bins(b) & input < bins(b+1));      
   out(ind) = level(b);
end
end