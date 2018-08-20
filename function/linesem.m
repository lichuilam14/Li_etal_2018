function [handies]= linesem(in,varargin)

colorcycle={[0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0      0.4470 0.7410],[0  0.4533 0.4533],...
    [0.8500 0.3250 0.0980],[0.2266 0.2266 1],[0.4940 0.1840 0.5560],[0.9066 0  0.4533],...
    [0 0 0]...
    [0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0      0.4470 0.7410],[0  0.4533 0.4533],...
    [0.8500 0.3250 0.0980],[0.2266 0.2266 1],[0.4940 0.1840 0.5560],[0.9066 0  0.4533],...
    [0 0 0]...
    [0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0      0.4470 0.7410],[0  0.4533 0.4533],...
    [0.8500 0.3250 0.0980],[0.2266 0.2266 1],[0.4940 0.1840 0.5560],[0.9066 0  0.4533],...
    [0 0 0]...
    [0.6350 0.0780 0.1840],[0.4660 0.6740 0.1880],[0      0.4470 0.7410],[0  0.4533 0.4533],...
    [0.8500 0.3250 0.0980],[0.2266 0.2266 1],[0.4940 0.1840 0.5560],[0.9066 0  0.4533],...
    [0 0 0]};


hold on;
linespecs =[];
barspecs  =[];
xspecs    =[];
alph      =[];
barwid    =[];
condnames =[];
CIfactor  =1;

% get variable arguments
if length(varargin)>0,linespecs=varargin{1};,end
if length(varargin)>1,barspecs=varargin{2};,end
if length(varargin)>2,xspecs=varargin{3};,end
if length(varargin)>3,condnames=varargin{4};,end
if length(varargin)>4,barwid=varargin{5};,end
if length(varargin)>5,alph=varargin{6};,end

if isempty(barwid),barwid=100;,end
if isempty(barspecs),barspecs={'color',[0.5 0.5 0.5],'linewidth',1};,end
if isempty(linespecs)
    linespecs={};
end

if isempty(linespecs)
    x=iscell(in);
    for n=1:size(in,2-x);
        linespecs{n}={'color',colorcycle{n},'linewidth',2,'LineStyle','-'};
    end
end



if ~iscell(in)
    
    if ndims(in)==2 | size(in,2)==1;
        for n=1:2;
            in1(:,n,:)=in;
        end
        in=in1;
    end
    
    stein=squeeze(nanstd(in,[],1))/sqrt(size(in,1));
    meanin=squeeze(nanmean(in));
    CIup  = meanin + (stein*CIfactor);
    CIlw  = meanin - (stein*CIfactor);
    
    np=size(in,3);
    nc=size(in,2);
    
else
    
    if ndims(in)==1;
        for n=1:length(in);
            meanin(n,:)=nanmean(in{n});
            stein(n,:)=squeeze(nanstd(in{n},[],1))/sqrt(size(in{n},1));
            CIup(n,:)  = meanin(n,:) + (stein(n,:)*CIfactor);
            CIlw(n,:)  = meanin(n,:) - (stein(n,:)*CIfactor);
            
        end
        
        meanin=[meanin';meanin'];
        stein=[stein';stein'];
        CIup =  [CIup';CIup'];
        CIlw =  [CIlw';CIlw'];
        
        nc=1;
        np=length(in);
        
        
    elseif ndims(in)==2;
        for n=1:size(in,1);
            for m=1:size(in,2);
                meanin(n,m)=nanmean(in{n,m});
                stein(n,m)=squeeze(nanstd(in{n,m},[],1))/sqrt(size(in{n,m},1));
                CIup(n,m)  = meanin(n,m) + (stein(n,m)*CIfactor);
                CIlw(n,m)  = meanin(n,m) - (stein(n,m)*CIfactor);
            end
        end
        
        np=size(in,2);
        nc=size(in,1);
        
    end
end



ytick=ranger(meanin);

if ~isempty(alph)
    if size(in,2)>1 & size(in,3)>1;
        [F1 F2 F12 p1 p2 p12]=pti(in);
        if p1<alph | p2<alph | p12<alph
            disp(['omnibus: p1=',num2str(p1),', p2=',num2str(p2),', p12=',num2str(p12)]);
            if nc==2; % if 2 conditiions, do t-test
                for n=1:size(in,3);
                    [ho pval(n) rci stats]=ttest_t(in(:,1,n)-in(:,2,n),0,alph);
                    if ho>0;
                        line([n,n],[ytick(1),ytick(end)],'linewidth',barwid,'color',[0.9 0.9 0.9]);
                    end
                end
                disp(['min p val:',num2str(min(pval))]);
            else  %if >2 conditiions, do anova
                for n=1:size(in,3);
                    [f1 f2 f12 p1 p2 p12]=pti(in(:,:,n));
                    if p1<alph;
                        line([n,n],[ytick(1),ytick(end)],'linewidth',barwid,'color',[0.9 0.9 0.9]);
                    end
                end
            end
        else
            disp(['n.s. p1=',num2str(p1),', p2=',num2str(p2),', p12=',num2str(p12)]);
        end
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% insert barspecs info %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toeval{1}='line([p,p],[CIlw(c,p),CIup(c,p)]';
toeval{2}='line([p-1/barwid,p+1/barwid],[CIlw(c,p),CIlw(c,p)]';
toeval{3}='line([p-1/barwid,p+1/barwid],[CIup(c,p),CIup(c,p)]';

for statement=1:length(toeval);
    for e=1:length(barspecs);
        if ischar(barspecs{e});
            toeval{statement}=[toeval{statement},',','''',barspecs{e},''''];
        elseif isnumeric(barspecs{e});
            if length(barspecs{e})>1;
                toeval{statement}=[toeval{statement},',[',num2str(barspecs{e}),']'];
            else
                toeval{statement}=[toeval{statement},',',num2str(barspecs{e})];
            end
        end
    end
    toeval{statement}=[toeval{statement},');'];
    %         disp(toeval{statement})    %debugging
end

longeval='for c=1:nc;for p=1:np;';
for statement=1:length(toeval);
    longeval=[longeval toeval{statement} ';'];
end

longeval=[longeval ';end;end'];

try
    eval(longeval);
catch
    disp('bars didnt work');
    disp(longeval);    %debugging
end

if ~isempty(linespecs);
    if ~iscell(linespecs{1});
        for c=1:nc;
            linespecs1{c}=linespecs;
        end
        linespecs=linespecs1;
    end
end



for c=1:nc;
    
    meanin_t=meanin';
    
    lineeval=['hog = plot(meanin_t(:,c)'];
    doteval =['plot(meanin_t(:,c),'];
    
    dot   = '.';
    clr   = 'color';
    mark  = 'markersize';
    DOT   = strcat('''',dot,'''');
    CLR   = strcat('''',clr,'''');
    MARK  = strcat('''',mark,'''');
    
    
    
    if ~isempty(linespecs);
        for e=1:length(linespecs{c});
            if ischar(linespecs{c}{e});
                lineeval=[lineeval,',','''',linespecs{c}{e},''''];
                doteval = [doteval,DOT,',',CLR,',[',num2str(linespecs{c}{2}),'],',MARK,',10'];
            elseif isnumeric(linespecs{c}{e})
                if length(linespecs{c}{e})>1;
                    lineeval=[lineeval,',[',num2str(linespecs{c}{e}),']'];
                    doteval = [doteval,DOT,',',CLR,',[',num2str(linespecs{c}{2}),'],',MARK,',10'];
                else
                    lineeval=[lineeval,',',num2str(linespecs{c}{e})];
                    doteval = [doteval,DOT,',',CLR,',[',num2str(linespecs{c}{2}),'],',MARK,',10'];
                end
            end
        end
    end
    
    lineeval=[lineeval,');'];
    doteval = [doteval,');'];
    leg = lineeval;
    L = doteval;
    try
        eval(lineeval);
        eval(doteval);
    catch
        disp('line or dot didnt work');
    end
    
    handies(c)=hog;
end

if ~isempty(condnames);
    legend(handies,condnames,'location','northeast');legend boxoff;
end

for statements=1:length(xspecs);
    if isnumeric(xspecs{statements}{2})
        xeval=['set(gca,''',xspecs{statements}{1},'''',',[','',num2str(xspecs{statements}{2}),']);'];
    else
        xeval=['set(gca,''',xspecs{statements}{1},''',{''',xspecs{statements}{2}{1},''',''',xspecs{statements}{2}{2},'''});'];
    end
    
    try
        eval(xeval);
    catch
    end
end




end