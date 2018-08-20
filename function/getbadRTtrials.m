function out = getbadRTtrials(sub,RT,q1,q2)
uSub = unique(sub);
out = [];
for s = 1:length(uSub)
    sRT = RT(sub == uSub(s));
    badRT = sRT>quantile(sRT,[q1])|sRT<quantile(sRT,[q2]);
    out  = [out badRT];
end
end