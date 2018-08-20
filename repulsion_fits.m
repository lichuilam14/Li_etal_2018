clear
clc
addpath('function');
addpath('jbtools');
allhRT = [];
allmRT = [];
for exp = 1:2
    switch exp
        case 1
            clearvars -except allhRT allmRT allbestmaxtw allbestSN
            load('data/PNAS_Exp2a.mat');
            fig_num = 'Fig. S4';
        case 2
            clearvars -except allhRT allmRT allbestmaxtw allbestSN
            load('data/PNAS_Exp2b.mat');
            fig_num = 'Fig. S5';
    end
    %% variables
    badRTtrials = getbadRTtrials(data.sub,data.RT,0.99,0.01);
    index       = ~badRTtrials & data.cor ==1;
    meanRT      = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),abs(data.targetMean(index)),abs(data.flankerMean(index)),data.congruency(index));
    uSub        = unique(data.sub);
    uFV         = unique(data.flankerSTD);
    uFM         = unique(abs(data.flankerMean));
    uTM         = unique(abs(data.targetMean));
    cong        = unique(data.congruency);
    numcond     = length(uFV)*length(uFM)*length(uTM)*length(cong);
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
    for s = 1:length(uSub)
        for v = 1:length(uFV)
            for f = 1:length(uFM)
                fm = uFM(f);
                for t = 1:length(uTM)
                    tm = uTM(t);
                    for c = 1:length(cong)
                        if c==2; fm=-fm;end
                        mR(v,t,f,c) = popcode(tm,fm,uFV(v)+bestSN(s),fspace,bestmaxtw(s));
                    end
                end
            end
        end
        R               = 1./reshape(abs(mR),[1 numcond]);
        hRT             = reshape(squeeze(meanRT(s,:,:,:,:)),[1 numcond]);
        [b dev]         = glmfit(abs(R),hRT);
        modelRT(s,:)    = b(1) + b(2).*abs(R);
    end
    %% plot (separate)
    % reshape
    reshape_mRT         = reshape(modelRT,[length(uSub) length(uFV) length(uTM) length(uFM) length(cong)]);
    rmRT(:,:,:,[3 2 1]) = reshape_mRT(:,:,:,:,1);
    rmRT(:,:,:,[4:6])   = reshape_mRT(:,:,:,:,2);
    rhRT(:,:,:,[3 2 1]) = meanRT(:,:,:,:,1);
    rhRT(:,:,:,[4:6])   = meanRT(:,:,:,:,2);
    allhRT              = [allhRT;rhRT];
    allmRT              = [allmRT;rmRT];
    % define color
    ci = linspace(0,1,length(uFM));
    for t = 1:length(uTM)
        for f = 1:length(uFM)*2
            if f <4;tcolmat{t,f} = [1,1-ci(f),0];end
            if f >3;tcolmat{t,f} = [0,1-ci(f-3),1];end
        end
    end
    
    % surface plot
    f = figure;
    set(f,'Units','inches','position',[0,0,10,10])
    for v = 1:3
        subplot(3,3,v);
        imagesc(squeeze(mean(rhRT(:,v,:,:),1)));
        caxis([0.5 0.62])
        colormap jet
        xticks([1:6])
        xticklabels({'H','M','L','L','M','H'})
        xlabel('Flanker mean orientation')
        line([3.5 3.5],[0 4],'color','k','LineStyle',':','LineWidth',3)
        yticks([1:3])
        yticklabels({'L','M','H'})
        ylabel('Target orientation')
        text(0.6,3.3,'Cong');
        text(5.2,3.3,'Incong');
        
        subplot(3,3,v+3);
        imagesc(squeeze(mean(rmRT(:,v,:,:),1)));
        caxis([0.5 0.62])
        colormap jet
        xticks([1:6])
        xticklabels({'H','M','L','L','M','H'})
        xlabel('Flanker mean orientation')
        line([3.5 3.5],[0 4],'color','k','LineStyle',':','LineWidth',3)
        yticks([1:3])
        yticklabels({'L','M','H'})
        ylabel('Target orientation')
        text(0.6,3.3,'Cong');
        text(5.2,3.3,'Incong');
    end
    
    % cross plot
    for v = 1:3
        subplot(3,3,v+6)
        hold on;
        line([0.52 0.62],[0.52,0.62],'color','k');
        for t = 1:3
            for f = 1:length(uFM)*2
                tcol = tcolmat{t,f};
                if t ==1;mark = 'o';end
                if t ==2;mark = 'v';end
                if t ==3;mark = 's';end
                hRT_m    = squeeze(mean(rhRT(:,v,t,f),1));
                hRT_ste  = ste(squeeze(rhRT(:,v,t,f)));
                hRT_CI   = [hRT_m-hRT_ste hRT_m+hRT_ste];
                mRT_m    = squeeze(mean(rmRT(:,v,t,f),1));
                mRT_ste  = ste(squeeze(rmRT(:,v,t,f)));
                mRT_CI   = [mRT_m-mRT_ste mRT_m+mRT_ste];
                hold on;
                line([hRT_CI],[mRT_m mRT_m],'Color',tcol,'linewidth',1.2);
                line([hRT_m hRT_m],[mRT_CI],'Color',tcol,'linewidth',1.2);
            end
        end
        
        for t = 1:3
            for f = 1:length(uFM)*2
                tcol = tcolmat{t,f};
                if t ==1;mark = 'o';end
                if t ==2;mark = 'v';end
                if t ==3;mark = 's';end
                hRT_m    = squeeze(mean(rhRT(:,v,t,f),1));
                mRT_m    = squeeze(mean(rmRT(:,v,t,f),1));
                plot(hRT_m,mRT_m,mark,'markerfacecolor',tcol,'markersize',10,'linewidth',1.2,'markeredgecolor','k');
            end
        end
        xlim([0.52 0.62]);ylim([0.52 0.62])
        xlabel('human RT(s)');ylabel('model RT(s)')        
    end
    suptitle(fig_num);
end

%% plot (comcatenated)
% surface plot
f = figure;
set(f,'Units','inches','position',[0,0,10,10])
for v = 1:3
    subplot(3,3,v);
    imagesc(squeeze(mean(allhRT(:,v,:,:),1)));
    caxis([0.5 0.62])
    colormap jet
    xticks([1:6])
    xticklabels({'H','M','L','L','M','H'})
    xlabel('Flanker mean orientation')
    line([3.5 3.5],[0 4],'color','k','LineStyle',':','LineWidth',3)
    yticks([1:3])
    yticklabels({'L','M','H'})
    ylabel('Target orientation')
    text(0.6,3.3,'Cong');
    text(5.2,3.3,'Incong');
    
    subplot(3,3,v+3);
    imagesc(squeeze(mean(allmRT(:,v,:,:),1)));
    caxis([0.5 0.62])
    colormap jet
    xticks([1:6])
    xticklabels({'H','M','L','L','M','H'})
    xlabel('Flanker mean orientation')
    line([3.5 3.5],[0 4],'color','k','LineStyle',':','LineWidth',3)
    yticks([1:3])
    yticklabels({'L','M','H'})
    ylabel('Target orientation')
    text(0.6,3.3,'Cong');
    text(5.2,3.3,'Incong');
end

% cross plot
for v = 1:3
    subplot(3,3,v+6)
    hold on;
    line([0.52 0.62],[0.52,0.62],'color','k');
    for t = 1:3
        for f = 1:length(uFM)*2
            tcol = tcolmat{t,f};
            if t ==1;mark = 'o';end
            if t ==2;mark = 'v';end
            if t ==3;mark = 's';end
            hRT_m    = squeeze(mean(allhRT(:,v,t,f),1));
            hRT_ste  = ste(squeeze(allhRT(:,v,t,f)));
            hRT_CI   = [hRT_m-hRT_ste hRT_m+hRT_ste];
            mRT_m    = squeeze(mean(allmRT(:,v,t,f),1));
            mRT_ste  = ste(squeeze(allmRT(:,v,t,f)));
            mRT_CI   = [mRT_m-mRT_ste mRT_m+mRT_ste];
            hold on;
            line([hRT_CI],[mRT_m mRT_m],'Color',tcol,'linewidth',1.2);
            line([hRT_m hRT_m],[mRT_CI],'Color',tcol,'linewidth',1.2);
        end
    end
    
    for t = 1:3
        for f = 1:length(uFM)*2
            tcol = tcolmat{t,f};
            if t ==1;mark = 'o';end
            if t ==2;mark = 'v';end
            if t ==3;mark = 's';end
            hRT_m    = squeeze(mean(allhRT(:,v,t,f),1));
            mRT_m    = squeeze(mean(allmRT(:,v,t,f),1));
            plot(hRT_m,mRT_m,mark,'markerfacecolor',tcol,'markersize',10,'linewidth',1.2,'markeredgecolor','k');
        end
    end
    xlim([0.52 0.62]);ylim([0.52 0.62])
    xlabel('human RT(s)');ylabel('model RT(s)')
end
suptitle('Fig. 6')