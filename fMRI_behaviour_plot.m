clear 
clc
load('data/PNAS_Exp3_fMRI.mat');
addpath('function');
addpath('jbtools');
%% variables
Tm           = data.targetMean;
Fm           = data.flankerMean;
Fv           = data.flankerSTD;
uSub         = unique(data.sub);
badRTtrials  = getbadRTtrials(data.sub,data.RT,0.99,0.01);
ES_maxtw     = linspace(5,60,30);
ES_SN        = linspace(0,10,20);
fspace       = -90:90;
%%  gain model
for s= 1:length(uSub)
    progressbar(length(uSub),s);
    clear indx subRT R
    indx = find(data.sub==uSub(s) & data.cor==1 &~badRTtrials);
    subRT = data.RT(indx);
    for mm = 1:length(ES_maxtw)
        maxtw = ES_maxtw(mm);
        for nn = 1:length(ES_SN)
            sn = ES_SN(nn);
            for t = 1:length(indx)             
                R(t) = popcode(Tm(indx(t)),Fm(indx(t)),Fv(indx(t))+sn,fspace,maxtw);
            end
            [b dev]             = glmfit(1./abs(R(:)),subRT(:));
            allparam(s,mm,nn,:) = b;            
            SSE(s,mm,nn)        = nansum(((b(1) + b(2)./abs(R(:))) - subRT(:)).^2);
        end
    end    
end
for s = 1:length(uSub);
    subSSE = squeeze(SSE(s,:,:));
    [bestmaxtw(s) bestSN(s)] = find(subSSE == min(subSSE(:)));
end
allmRT     = [];
for s = 1:length(uSub)
    clear gainR mRT 
    indx     = find(data.sub == uSub(s));
    p        = allparam(s,bestmaxtw(s),bestSN(s),:);
    for t = 1:length(indx)
        R              = popcode(Tm(indx(t)),Fm(indx(t)),Fv(indx(t))+ ES_SN(bestSN(s)),fspace,ES_maxtw(bestmaxtw(s)));
        R(abs(R)< 0.1) = 0.1.*sign(R(abs(R)<0.1));
        mRT(t)         = p(1) + p(2)./abs(R);
    end
    allmRT = [allmRT mRT]; 
end
%% plot 
Tm      = abs(data.targetMean);
Fm      = abs(data.flankerMean);
Fv      = data.flankerSTD;
Tmlevel = discreter(Tm,15,30);
Fmlevel = discreter(Fm,15,30);
Fvlevel = discreter(Fv,15);
uTm     = unique(Tmlevel);
uFm     = unique(Fmlevel);
uFv     = unique(Fvlevel);
cong    = unique(data.congruency);
index   = data.cor == 1 & ~badRTtrials;
clear mRT hRT
for s = 1:length(uSub)
    for v = 1:length(uFv)
        for t = 1:length(uTm)
            for f = 1:length(uFm)
                for c = 1:length(cong)
                    ind = find(data.sub== uSub(s)& data.congruency==cong(c) &Tmlevel == uTm(t)& Fmlevel == uFm(f)& Fvlevel==uFv(v) & index);
                    hRT(s,v,t,f,c) = nanmean(data.RT(ind));
                    mRT(s,v,t,f,c)= nanmean(allmRT(ind));
                end
            end
        end
    end
end
rhRT(:,:,:,[3 2 1]) = hRT(:,:,:,:,1);
rhRT(:,:,:,[4:6])   = hRT(:,:,:,:,2);
rmRT(:,:,:,[3 2 1]) = mRT(:,:,:,:,1);
rmRT(:,:,:,[4:6])   = mRT(:,:,:,:,2);
ci = linspace(0,1,length(unique(uFm)));
for t = 1:3
    for f = 1:length(unique(uFm))*2
        if f <4;tcolmat{t,f} = [1,1-ci(f),0];end
        if f >3;tcolmat{t,f} = [0,1-ci(f-3),1];end
    end
end

fig = figure;
set(fig,'Units','inches','position',[0,0,7,3])
suptitle('Fig.S8')
for fv = 1:2
    fig_n = {'Low Flanker std','Low Flanker std'};
    subplot(1,2,fv)
    line([0.52 0.8],[0.58 0.7],'Color','k');
    xlim([0.52 0.8]);ylim([0.58 0.7]);
    xlabel('Human RT (s)')
    ylabel('Model RT(s)');
    
    for t = 1:3
        for f = 1:length(unique(uFm))*2
            tcol = tcolmat{t,f};
            if t ==1;mark = 'o';end
            if t ==2;mark = 'v';end
            if t ==3;mark = 's';end
            hRT_m    = squeeze(nanmean(rhRT(:,fv,t,f),1));
            hRT_ste  = nanste(nanmean(rhRT(:,fv,t,f)));
            hRT_CI   = [hRT_m-hRT_ste hRT_m+hRT_ste];
            mRT_m    = squeeze(nanmean(rmRT(:,fv,t,f),1));
            mRT_ste  = nanste(squeeze(rmRT(:,fv,t,f)));
            mRT_CI   = [mRT_m-mRT_ste mRT_m+mRT_ste];
            hold on;
            line([hRT_CI],[mRT_m mRT_m],'Color',tcol,'linewidth',1.2);
            line([hRT_m hRT_m],[mRT_CI],'Color',tcol,'linewidth',1.2);
        end
    end
    for t = 1:3
        for f = 1:length(unique(uFm))*2
            tcol = tcolmat{t,f};
            if t ==1;mark = 'o';end
            if t ==2;mark = 'v';end
            if t ==3;mark = 's';end
            hRT_m    = squeeze(nanmean(rhRT(:,fv,t,f),1));
            mRT_m    = squeeze(nanmean(rmRT(:,fv,t,f),1));
            hold on;
            plot(hRT_m,mRT_m,mark,'markerfacecolor',tcol,'markersize',10,'linewidth',1.2,'markeredgecolor','k');
        end
    end      
end