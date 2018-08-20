clear
clc
addpath('function')
%% conflict model 2
allmeanRT = [];allhmeanRT = [];
for exp = 1:2
    switch exp
        case 1          
            clearvars -except allmeanRT allhmeanRT
            load('data/PNAS_Exp2a.mat')
            tm  = data.targetMean;
            fm  = data.allangles(2:end,:);                        
        case 2
            clearvars -except allmeanRT allhmeanRT
            load('data/PNAS_Exp2b.mat')
            tm  = data.targetMean.*180;
            fm  = data.allcolor(2:end,:).*180;
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
        indx   = find(data.sub==uSub(s) & index);
        subtm  = tm(indx);
        suballf= fm(:,indx); 
        subRT  = data.RT(indx);
        for mm = 1:length(ES_weight)
            w  = ES_weight(mm);
            for t = 1:length(indx)
                currt        = subtm(t);
                currf        = suballf(:,t);                
                eviL         = (abs(currt).*(currt<0).*w)  + sum((abs(currf).*(currf<0).*((1-w)./6)));
                eviR         = (abs(currt).*(currt>0).*w)  + sum((abs(currf).*(currf>0).*((1-w)./6)));
                conflict(t)  = eviL.*eviR;
            end
            
            conf     = abs(conflict(:));
            bc       = glmfit(conf,subRT(:));
            RT_c     = bc(1)+conf*bc(2);
            SSE_c(s,mm)  = sum((RT_c - subRT').^2);            
        end
        clear indx subtm suballf subRT RT_c conf conflict
        progressbar(length(uSub),s);
    end
    %% find best
    for s = 1:length(uSub)
        bw = find(SSE_c(s,:)==min(SSE_c(s,:)));
        bestweight(s)= bw(1);
    end
    %% best parameterisation
    allmRT = [];
    for s = 1:length(uSub)
        indx   = find(data.sub==uSub(s));
        subtm  = tm(indx);
        suballf= fm(:,indx);
        subRT  = data.RT(indx);
        w     = ES_weight(bestweight(s));
        for t = 1:length(indx)
            currt        = subtm(t);
            currf        = suballf(:,t);
            eviL         = (abs(currt).*(currt<0).*w)  + sum((abs(currf).*(currf<0).*((1-w)./6)));
            eviR         = (abs(currt).*(currt>0).*w)  + sum((abs(currf).*(currf>0).*((1-w)./6)));
            conflict(t)  = eviL.*eviR;
        end                
        conf = abs(conflict(:));
        bc = glmfit(conf,subRT(:));
        mRT = bc(1) + bc(2).*conf;
        allmRT = [allmRT;mRT];
    end
    %% mean model RT
    index     = ~eventrials & ~badRTtrials & data.cor ==1;    
    mean_mRT  = jb_getvector(allmRT(index),data.sub(index),data.flankerSTD(index),abs(data.targetMean(index)),abs(data.flankerMean(index)),data.congruency(index));
    mean_hRT  = jb_getvector(data.RT(index),data.sub(index),data.flankerSTD(index),abs(data.targetMean(index)),abs(data.flankerMean(index)),data.congruency(index));
    allmeanRT = [allmeanRT;mean_mRT];
    allhmeanRT= [allhmeanRT;mean_hRT];
end

%% plot
rhRT(:,:,:,[3 2 1]) = allhmeanRT(:,:,:,:,1);
rhRT(:,:,:,[4:6])   = allhmeanRT(:,:,:,:,2);
rmRT(:,:,:,[3 2 1]) = allmeanRT(:,:,:,:,1);
rmRT(:,:,:,[4:6])   = allmeanRT(:,:,:,:,2);
ci = linspace(0,1,size(rmRT,4)/2);
for t = 1:3
    for f = 1:size(rmRT,4)
        if f <4;tcolmat{t,f} = [1,1-ci(f),0];end
        if f >3;tcolmat{t,f} = [0,1-ci(f-3),1];end
    end
end
titlemat  ={'Low','Medium','High'};
f = figure;
set(f,'Units','inches','position',[0,0,12,3])
for v = 1:3
    subplot(1,3,v)
    hold on;
    line([0.52 0.62],[0.52,0.62],'color','k');
    for t = 1:3
        for f = 1:6
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
        for f = 1:6
            tcol = tcolmat{t,f};
            if t ==1;mark = 'o';end
            if t ==2;mark = 'v';end
            if t ==3;mark = 's';end
            hRT_m    = squeeze(mean(rhRT(:,v,t,f),1));
            mRT_m    = squeeze(mean(rmRT(:,v,t,f),1));
            plot(hRT_m,mRT_m,mark,'markerfacecolor',tcol,'markersize',10,'linewidth',1.2,'markeredgecolor','k');
        end
    end
    xlim([0.52 0.62]);ylim([0.52 0.62]);
    xlabel('human RT(s)');ylabel('model RT(s)');
    title([titlemat{v},'flanker variance']);
end




