function [wwX, wherex] =barsem(x,varargin);
wwX      =[];
plottype ='b';
inwid    =20;
stecolor1='k';
stecolor2='w';
steshape ='o';
steline1 =1;
steline2 =1;
mks      =7;

if length(varargin)>0,plottype=varargin{1};,end
if length(varargin)>1,inwid=varargin{2};,end

%% get colours
if length(varargin)>2,
    stecolors=varargin{3};
    if length(stecolors)==2;
        stecolor1=stecolors{1};
        stecolor2=stecolors{2};
    else
        stecolor1=stecolors;
        stecolor2=stecolors;
    end
end
if length(varargin)>3,steshape=varargin{4};,end

if length(varargin)>4,
    stelines=varargin{5};
    if length(stecolors)==2;
        steline1=stelines{1};
        steline2=stelines{2};
    else
        steline1=stelines;
        steline2=stelines;
    end
end



% ci =linspace(0,1,10);
% for i = 1:length(ci)
%     col(i,:) = [1-ci(i) 0 0+ci(i)];
% end

col(1,:) = [80,80,80]./255;
col(2,:) = [255,255,255]./255;
col(3,:) = [200,200,200]./255;

colormap(col)
% colormap colorcube
if ndims(x)>3
    error('cant plot a 4d matrix');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FIRST DEAL WITH 2-D DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ndims(x)==2;
    
    xlim([0 size(x,2)+1]);
    
    np=size(x,2);
    
    % get standard errors
    steX=(nanstd(x,[],1))/sqrt(size(x,1));     %standard error is std/root n
    
    % get mean
    mx=squeeze(nanmean(x));
    
    
    
    
    hold on;
    
    
    h = bar(mx);
    set(gca,'clim',[0 1]);
    wwX = h;
    wherex = bsxfun(@plus , h.XData,[h.XOffset]);
    
    
    hold on
    %plot error bars first
    for n=1:np;
        line([n,n],[mx(n)-(steX(n).*1.96),mx(n)+(steX(n).*1.96)],'color',[0.5 0.5 0.5],'linewidth',steline1);
        if size(mx,2)<100;
            line([n-1/inwid,n+1/inwid],[mx(n)-(steX(n).*1.96),mx(n)-(steX(n).*1.96)],'color',[0.5 0.5 0.5],'linewidth',steline1);
            line([n-1/inwid,n+1/inwid],[mx(n)+(steX(n).*1.96),mx(n)+(steX(n).*1.96)],'color',[0.5 0.5 0.5],'linewidth',steline1);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% WASN'T THAT EASY? WHAT ABOUT 3-D DATA? %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif ndims(x)==3;
    xlim([0 size(x,3)+1]);
        
    bw=1;
    mx=squeeze(nanmean(x));
    
    h=bar(mx');
    h1=get(gca,'Children');a=get(h1);;
    np=size(x,3);
    X=x(:,:);
    steX=nanstd(X,[],1)/sqrt(size(x,1));
    CIX = steX * 1.96;
    steX=reshape(steX,size(x,2),np);
    CIX =reshape(CIX,size(x,2),np);
    wwX = h;
    for l=1:size(x,2);
        wherex(l,:) = bsxfun(@plus , h(l).XData,[h(l).XOffset]);
        wherey(l,:) = h(l).YData;
        wid = 0.01;
        for n=1:np;
            line([wherex(l,n) wherex(l,n)],[wherey(l,n)-CIX(l,n) wherey(l,n)+CIX(l,n)],'color',[0.5 0.5 0.5],'linewidth',bw);
            line([wherex(l,n)-wid wherex(l,n)+wid],[wherey(l,n)-CIX(l,n) wherey(l,n)-CIX(l,n)],'color',[0.5 0.5 0.5],'linewidth',bw);
            line([wherex(l,n)-wid wherex(l,n)+wid],[wherey(l,n)+CIX(l,n) wherey(l,n)+CIX(l,n)],'color',[0.5 0.5 0.5],'linewidth',bw);
        end
    end
    
    
    
end
