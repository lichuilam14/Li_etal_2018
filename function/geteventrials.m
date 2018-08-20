function out = geteventrials(condition)
uCond  =  unique(condition);
allind = [];
for c = 1:length(uCond)
    ind        = find(condition==uCond(c));
    allind     = [allind ind];
end
out(allind(2:2:end)) = 1;
end