clear
clc
addpath('function')
%% conflict model 1
f = figure;
set(f,'Units','inches','position',[0,0,7,3])
expt = {'Exp. 1a','Exp. 1b'};
for exp = 1:2
    switch exp
        case 1
            load('data/PNAS_Exp1a.mat')
        case 2
            load('data/PNAS_Exp1b.mat')
    end
    %% get mean RT
    eventrials  = geteventrials(data.flankerSTD);
    badRTtrials = getbadRTtrials(data.sub,data.RT,0.99,0.01);
    uSub        = unique(data.sub);
    uFV         = unique(data.flankerSTD);
    cong        = unique(data.congruency);
    index       = eventrials & ~badRTtrials & data.cor ==1;
    meanRT      = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),data.congruency(index));
    %% parameters
    ES_weight   = linspace(0,1,50);
    %% model
    for s = 1:length(uSub)
        subRT = squeeze(meanRT(s,:,:));
        for mm = 1:length(ES_weight)
            w  = ES_weight(mm);
            for fv = 1:length(uFV)
                for c = 1:length(cong)
                    if c == 1; fm = 45 ;end
                    if c == 2; fm = -45;end
                    eviL           = abs(fm).*(fm<0).*w;
                    eviR           =  (45) .*(1-w) + (abs(fm).*(fm>0).*w);
                    conflict(fv,c) = eviL.*eviR;
                end
            end
            conf = abs(conflict(:));
            bc = glmfit(conf,subRT(:));
            RT_c = bc(1)+conf*bc(2);
            SSE_c(s,mm) = sum((RT_c-subRT(:)).^2);
        end
        progressbar(length(uSub),s);
    end
    %% find best
    for s = 1:length(uSub)
        bw = find(SSE_c(s,:)==min(SSE_c(s,:)));
        bestweight(s)= bw(1);
    end
    %% best parameterisation
    index     = ~eventrials & ~badRTtrials & data.cor ==1;
    meanRT    = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),data.congruency(index));
    RT_c        = nan(length(uSub),3,2);
    for s = 1:length(uSub)
        subRT = squeeze(meanRT(s,:,:));
        w     = ES_weight(bestweight(s));
        for fv = 1:length(uFV)
            for c = 1:length(cong)
                if c == 1; fm = 45 ;end
                if c == 2; fm = -45;end
                eviL           = abs(fm).*(fm<0).*w;
                eviR           =  (45) .*(1-w) + (abs(fm).*(fm>0).*w);
                conflict(fv,c) = eviL.*eviR;
            end
        end
        conf = abs(conflict(:));
        bc = glmfit(conf,subRT(:));
        RT_c(s,:,:) = reshape(bc(1)+conf*bc(2),[length(uFV) length(cong)]);
    end
    %% plot
    subplot(1,2,exp);
    linesem(RT_c.*1000);
    title([expt{exp}]);
    xticks([1 2])
    xticklabels({'Cong','Incong'})
    ylabel('RT(ms)')
    clc;
end
suptitle('Conflict model 1');


%% conflict model 2
clear
clc
f = figure;
set(f,'Units','inches','position',[0,0,7,3])
expt = {'Exp. 1a','Exp. 1b'};
for exp = 1:2
    switch exp
        case 1
            load('data/PNAS_Exp1a.mat')
        case 2
            load('data/PNAS_Exp1b.mat')
    end
    %% get mean RT
    eventrials  = geteventrials(data.flankerSTD);
    badRTtrials = getbadRTtrials(data.sub,data.RT,0.99,0.01);
    uSub        = unique(data.sub);
    index       = eventrials & ~badRTtrials & data.cor ==1;
    %% parameters
    ES_weight   = linspace(0,1,50);
    %% model
    for s = 1:length(uSub)
        indx  = find(data.sub == uSub(s) & index);
        subRT = data.RT(indx);
        subtm = data.targetMean(indx);
        subf  = data.allangles(2:end,indx);
        for mm = 1:length(ES_weight)
            w  = ES_weight(mm);
            eviL           = (abs(subtm).*(subtm<0).*w)  + sum((abs(subf).*(subf<0).*((1-w)./6)));
            eviR           = (abs(subtm).*(subtm>0).*w)  + sum((abs(subf).*(subf>0).*((1-w)./6)));
            conflict        = eviL.*eviR;
            conf = abs(conflict(:));
            bc = glmfit(conf,subRT(:));
            RT_c = bc(1)+conf*bc(2);
            SSE_c(s,mm) = sum((RT_c-subRT(:)).^2);
        end
        progressbar(length(uSub),s);
    end
    %% find best
    for s = 1:length(uSub)
        bw = find(SSE_c(s,:)==min(SSE_c(s,:)));
        bestweight(s)= bw(1);
    end
    %% best parameterisation
    allmRT    = [];
    for s = 1:length(uSub)
        indx  = find(data.sub == uSub(s));
        subRT = data.RT(indx);
        subtm = data.targetMean(indx);
        subf  = data.allangles(2:end,indx);
        w     = ES_weight(bestweight(s));
        eviL  = (abs(subtm).*(subtm<0).*w)  + sum((abs(subf).*(subf<0).*((1-w)./6)));
        eviR  = (abs(subtm).*(subtm>0).*w)  + sum((abs(subf).*(subf>0).*((1-w)./6)));
        conflict  = eviL.*eviR;
        conf      = abs(conflict(:));
        bc        = glmfit(conf,subRT(:));
        mRT       = bc(1)+ bc(2)*conf;
        allmRT    = [allmRT;mRT];
    end
    %% plot
    index     = ~eventrials & ~badRTtrials & data.cor ==1;
    RT_c      = jb_getvector(allmRT(index),data.sub(index),data.flankerSTD(index),data.congruency(index)).*1000;
    subplot(1,2,exp);
    linesem(RT_c.*1000);
    title([expt{exp}]);
    xticks([1 2])
    xticklabels({'Cong','Incong'})
    ylabel('RT(ms)')
    clc;
end
suptitle('Conflict model 2');
