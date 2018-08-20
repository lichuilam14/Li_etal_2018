clear
clc
addpath('function');
addpath('jbtools');
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
    fspace      = -90:90;
    ES_maxtw    = linspace(5,60,50);  
    ES_SN       = linspace(0,40,40);
    %% model
    for s = 1:length(uSub)        
        subRT = squeeze(meanRT(s,:,:));
        for mm = 1:length(ES_maxtw)
            maxtw = ES_maxtw(mm);
            for nn = 1:length(ES_SN)
                sn = ES_SN(nn);
                for fv = 1:3
                    for c = 1:2
                        if c == 1; fm = 45 ;end
                        if c == 2; fm = -45;end
                        mR(fv,c)       = popcode(45,fm,uFV(fv)+sn,fspace,maxtw);
                    end
                end
                R  = 1./abs(mR(:));
                b  = glmfit(R,subRT(:));
                RT_g = b(1)+R*b(2);
                SSE_g(s,mm,nn) = sum(abs(RT_g-subRT(:)));
            end
        end
        progressbar(length(uSub),s);
    end
    %% find best
    for s = 1:length(uSub)
        subSSE_g = squeeze(SSE_g(s,:,:));
        [bestmaxtw(s) bestsn(s)] = ind2sub(size(subSSE_g),find(subSSE_g==min(subSSE_g(:))));
    end
    %% best parameterisation
    index     = ~eventrials & ~badRTtrials & data.cor ==1;
    meanRT    = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),data.congruency(index));
    RT_g      = nan(length(uSub),3,2);
    for s = 1:length(uSub)
        subRT = squeeze(meanRT(s,:,:));
        maxtw = ES_maxtw(bestmaxtw(s));
        sn    = ES_SN(bestsn(s));
        for fv = 1:length(uFV)
            for c = 1:length(cong)
                if c == 1; fm = 45 ;end
                if c == 2; fm = -45;end
                mR(fv,c)       = popcode(45,fm,uFV(fv)+sn,fspace,maxtw);
            end
        end
        R  = 1./abs(mR(:));
        b  = glmfit(R,subRT(:));
        RT_g(s,:,:) = reshape(b(1)+R*b(2),[3 2]).*1000;
    end
    %% plot
    index     = data.cor ==1 & ~badRTtrials & rem(data.trial,2)==1;
    meanRT    = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),data.congruency(index));
    f = figure;
    set(f,'Units','inches','position',[0,0,3,3])
    [hh]= linesem(meanRT.*1000);
    hold on;plot(squeeze(mean(permute(RT_g(:,1,:),[1 3 2]),1)),'o','markersize',6,'linewidth',2.5,'markeredgecolor',[0.6350 0.0780 0.1840]);
    hold on;plot(squeeze(mean(permute(RT_g(:,2,:),[1 3 2]),1)),'o','markersize',6,'linewidth',2.5,'markeredgecolor',[0.4660 0.6740 0.1880]);
    hold on;plot(squeeze(mean(permute(RT_g(:,3,:),[1 3 2]),1)),'o','markersize',6,'linewidth',2.5,'markeredgecolor',[0      0.4470 0.7410]);
    legend([hh],{'low','medium','high'},'Location','southeast')
    xticks([1 2])
    xticklabels({'Cong','Incong'})
    ylabel('RT(ms)')
end






