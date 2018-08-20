clear
clc
addpath('function')
warning off
%%
allB1 = [];allB2 = [];allB3 = [];allB4 = [];
for exp = 1:2
    switch exp
        case 1
            clearvars -except allB*
            load('data/PNAS_Exp2a.mat');
            flankerMean = data.flankerMean;
            targetMean  = data.targetMean;
            flankerSTD  = data.flankerSTD;
            allflankers = data.allangles(2:end,:);            
            fig_num = 'Fig. S4';
        case 2
            clearvars -except allB*
            load('data/PNAS_Exp2b.mat');
            flankerMean = data.flankerMean.*180;
            targetMean  = data.targetMean.*180;            
            allflankers = data.allcolor(2:end,:).*180;            
            flankerSTD  = data.flankerSTD;
            fig_num = 'Fig. S5';
    end
    %% variables
    badRTtrials = getbadRTtrials(data.sub,data.RT,0.99,0.01);
    index       = ~badRTtrials & data.cor ==1;
    meanRT      = jb_getvector(data.RT(index),data.sub(index),flankerSTD(index),abs(targetMean(index)),abs(flankerMean(index)),data.congruency(index));
    uSub        = unique(data.sub);
    uFV         = unique(flankerSTD);
    uFM         = unique(abs(flankerMean));
    uTM         = unique(abs(targetMean));
    cong        = unique(data.congruency);
    numcond     = length(uFV)*length(uFM)*length(uTM)*length(cong);
    alltfdiff   = allflankers - targetMean;
    binspace    = [min(alltfdiff(:)) linspace(-60,60,10) max(alltfdiff(:))];
    %% parameters
    ES_maxtw    = linspace(5,100,60);
    ES_SN       = linspace(0.05,100,60);
    fspace      = linspace(-90,90,180);
    %% model
    for s = 1:length(uSub)
        for mm = 1:length(ES_maxtw)
            maxtw = ES_maxtw(mm);
            for nn = 1:length(ES_SN)
                sn = ES_SN(nn);
                for v = 1:length(uFV)
                    for f = 1:length(uFM)
                        fm = uFM(f);
                        for t = 1:length(uTM)
                            tm = uTM(t);
                            for c = 1:length(cong)
                                if c ==2;fm = -fm;end
                                mR(v,t,f,c) = popcode(tm,fm,uFV(v)+sn,fspace,maxtw);
                            end
                        end
                    end
                end
                R             = 1./reshape(abs(mR),[1 numcond]);
                hRT           = reshape(squeeze(meanRT(s,:,:,:,:)),[1 numcond]);
                [b dev]       = glmfit(abs(R),hRT);
                tempmRT       = b(1) + b(2).*abs(R);
                SSE(s,mm,nn)  = sum((tempmRT - hRT).^2);
            end
        end
        progressbar(length(uSub),s);
    end
    for s = 1:length(uSub)
        hRT          = reshape(meanRT(s,:,:,:,:),[1 numcond]);
        subSSE       = squeeze(SSE(s,:,:));
        [M N]        = ind2sub(size(subSSE),find(subSSE== min(subSSE(:))));
        bestmaxtw(s) = ES_maxtw(M);
        bestSN(s)    = ES_SN(N);
    end
    clear mR
    %% regression
    for s = 1:length(uSub)
        ind   = find(data.sub==uSub(s)& data.key_cat<2);
        maxtw = bestmaxtw(s);
        sn    = bestSN(s);
        subtm = targetMean(ind);
        subfm = flankerMean(ind);
        subfv = flankerSTD(ind);
        subRT = data.RT(ind);
        subkey= data.key_cat(ind);
        tfdiff= alltfdiff(:,ind);
        for t = 1:length(ind)
            mR(t)       = popcode(subtm(t),subfm(t),subfv(t)+sn,fspace,maxtw);
            mRT(t)      = 1./abs(mR(t));
            hkey(t)     = subkey(t);
            hRT(t)      = subRT(t);
            for n = 1:length(binspace)-1
                indx = find(tfdiff(:,t) >binspace(n) & tfdiff(:,t) <=binspace(n+1));
                XX(t,n) = sum(tfdiff(indx,t),1);
            end
        end
        mkey                       = sigmoidv(mR,0,1,0,5);
        pred(:,1)                  = subtm;
        pred(:,2:length(binspace)) = XX;
        B1(s,:)                    = glmfit(pred,hkey'==1,'binomial','link','probit');
        B2(s,:)                    = glmfit(pred,hRT');
        B3(s,:)                    = glmfit(pred,mkey'>0.5,'binomial','link','probit');
        B4(s,:)                    = glmfit(pred,mRT');
        
        clear h* m* tfdiff ind* pred badresp sub* XX
    end
    allB1 = [allB1;B1];
    allB2 = [allB2;B2];
    allB3 = [allB3;B3];
    allB4 = [allB4;B4];
end
%% plot
f = figure;
set(f,'Units','inches','position',[0,0,11,3])
subplot(1,4,1);
shadedErrorBar(1:length(binspace)-1,squeeze(mean(allB1(:,3:end),1)),squeeze(ste(allB1(:,3:end),1)),{'color',[0 0.8 0.8],'linewidth',3},1);xlim([1 11]);
subplot(1,4,2);
shadedErrorBar(1:length(binspace)-1,squeeze(mean(allB3(:,3:end),1)),squeeze(ste(allB3(:,3:end),1)),{'LineStyle',':','color',[0 0.8 0.8],'linewidth',3},1);xlim([1 11]);
subplot(1,4,3);
shadedErrorBar(1:length(binspace)-1,squeeze(mean(allB2(:,3:end),1)),squeeze(ste(allB2(:,3:end),1)),{'color',[0.6 0.4 1],'linewidth',3},1);xlim([1 11]);
subplot(1,4,4);
shadedErrorBar(1:length(binspace)-1,squeeze(mean(allB4(:,3:end),1)),squeeze(ste(allB4(:,3:end),1)),{'LineStyle',':','color',[0.6 0.4 1],'linewidth',3},1);xlim([1 11]);
warning on;
clc