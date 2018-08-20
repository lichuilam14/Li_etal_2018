clear
clc
addpath('function')
addpath('jbtools')
load('output/allgainR.mat');
load('output/allconflict.mat');
load('data/PNAS_Exp3_fMRI.mat');
model    = struct;
Tmlevel  = discreter(data.targetMean,-22.5,0,22.5);
Fmlevel  = discreter(data.flankerMean,-22.5,0,22.5);
Tm       = data.targetMean;
Tm(abs(Tm)< 0.1) = 0.1.*sign(Tm(abs(Tm)<0.1));
uSub     = unique(data.sub);
uTm      = unique(Tmlevel);
uFm      = unique(Fmlevel);
signallgainR = sign(allgainR);
absallgainR = abs(allgainR);
absallgainR(absallgainR<0.1) = 0.1;
allgainR = signallgainR.*absallgainR;
for s = 1:length(uSub)
    for t = 1:length(uTm)
        for f = 1:length(uFm)
            ind = find(data.sub==uSub(s)& Tmlevel == uTm(t)& Fmlevel ==uFm(f));
            model.R(s,t,f)        = squeeze(nanmean(1./abs(allgainR(ind))));
            model.Tm(s,t,f)       = squeeze(nanmean(abs(1./Tm(ind))));
            model.conflict(s,t,f) = squeeze(mean(abs(allconflict(ind))));
        end
    end
end
%%  plot
f   = figure;
set(f,'Units','inches','position',[0,0,10,3]);
fn      = fieldnames(model);
title_n = {'1/|R|','1/|Xi|','conflict'}

for n = 1:length(fn)
    subplot(1,length(fn),n)
    eval(['[hh] = linesem(model.',fn{n},')']);
    xlim([0.5 4.5]);
    xticks([1:4])
    xticklabels({'-High','-Low','+Low','+High'});
    xlabel('Flanker mean orientation')    
    title(title_n{n})
end
legend([hh],{'-High','-Low','+Low','+High'},'Location','best')