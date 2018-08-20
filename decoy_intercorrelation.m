%% fig 3 e-h intercorrelation of decoy effects
clear
clc
fspace = -90:90;
Xsmat  = 30;
X = [15,10];
Y = [10,15];
Z = [8.5,16.5];
maxtwmat = linspace(8,16,40);
subslope = ones(40,1).*0.08;
trials = 100;
iterations = 100;

fig =figure;
set(fig,'Units','inches','position',[0,0,80,3.5])
subplot(1,4,1)
line([0 25],[25 0],'Color',[0,0,0])
hold on;
plot(X(1),X(2),'ro','linewidth',5);
plot(Y(1),Y(2),'o','color',[0.5,0.5,0.5],'linewidth',5);
plot(Z(:,1),Z(:,2),'.','MarkerEdgeColor',[0,0.8,0.2],'markersize',30);

for n = 1:length(maxtwmat)
    maxtw = maxtwmat(n);
    scaly = 440;
    
    sigmat1 = maxtw-(normpdf(fspace,Z(1),Xsmat)*scaly);
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    
    sigmat2 = maxtw-(normpdf(fspace,Z(2),Xsmat)*scaly);
    PX2 = normpdf(X(2),fspace,sigmat2);
    PY2 = normpdf(Y(2),fspace,sigmat2);
    X2  = PX2*fspace';
    Y2  = PY2*fspace';
    allX = X1+X2;
    allY = Y1+Y2;
    XYdiff = allX-allY;
    
    for iter = 1:iterations
        ndiff = repmat(XYdiff,[trials,1])+randn(trials,1).*subslope(n);
        mch   =  ndiff>0;
        PRST_sim(n,iter) = mean(mch);
    end

end
% compromise
fspace = -90:90;
Xsmat  = 30;
X = [15,10];
Y = [10,15];
Z = [20,5];
plot(Z(:,1),Z(:,2),'.','MarkerEdgeColor',[0,1,1],'markersize',30);
for n = 1:length(maxtwmat)
    maxtw = maxtwmat(n);
    scaly = 440;
    
    sigmat1 = maxtw-(normpdf(fspace,Z(1),Xsmat)*scaly);
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    
    sigmat2 = maxtw-(normpdf(fspace,Z(2),Xsmat)*scaly);
    PX2 = normpdf(X(2),fspace,sigmat2);
    PY2 = normpdf(Y(2),fspace,sigmat2);
    X2  = PX2*fspace';
    Y2  = PY2*fspace';
    allX = X1+X2;
    allY = Y1+Y2;
    XYdiff = allX-allY;
    
    for iter = 1:iterations
        ndiff = repmat(XYdiff,[trials,1])+randn(trials,1).*subslope(n);
        mch   =  ndiff>0;
        PRST_comp(n,iter) = mean(mch);
    end

end

% attraction
fspace = -90:90;
Xsmat  = 30;
X = [15,10];
Y = [10,15];
a = sum(X);
Z = [13.5,8.5];

plot(Z(:,1),Z(:,2),'.','MarkerEdgeColor',[1,0,1],'markersize',30);
for n = 1:length(maxtwmat)
    maxtw = maxtwmat(n);
    scaly = 440;
    
    sigmat1 = maxtw-(normpdf(fspace,Z(1),Xsmat)*scaly);
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    
    sigmat2 = maxtw-(normpdf(fspace,Z(2),Xsmat)*scaly);
    PX2 = normpdf(X(2),fspace,sigmat2);
    PY2 = normpdf(Y(2),fspace,sigmat2);
    X2  = PX2*fspace';
    Y2  = PY2*fspace';
    allX = X1+X2;
    allY = Y1+Y2;
    XYdiff = allX-allY;
    for iter = 1:iterations
        ndiff = repmat(XYdiff,[trials,1])+randn(trials,1).*subslope(n);
        mch   =  ndiff>0;
        PRST_attr(n,iter) = mean(mch);
    end
end
ylim([0 25])
xlim([0 25])
xlabel('Attribute 1')
ylabel('Attribute 2')
%% plot
subplot(1,4,4)
suptitle('decoy intercorrelation - fig. 3')
for n = 1:length(maxtwmat)
    hold on;
    plot(PRST_sim(n,:),PRST_comp(n,:),'.','color',[1-(0.025*n),0,0.025*n],'markersize',5)
end

p2 = polyfit(PRST_sim(:),PRST_comp(:),1);
f2 = polyval(p2,PRST_sim(:));
hold on;plot(PRST_sim(:),f2,'color','k','linewidth',2);
xlabel('Similarity')
ylabel('Compromise')
ylim([0 1])
xlim([0 1])


subplot(1,4,2)
for n = 1:length(maxtwmat)
    hold on;
    plot(PRST_attr(n,:),PRST_comp(n,:),'.','color',[1-(0.025*n),0,0.025*n],'markersize',5)
end
xlabel('Attraction')
ylabel('Compromise')
p2 = polyfit(PRST_attr(:),PRST_comp(:),1);
f2 = polyval(p2,PRST_attr(:));
hold on;plot(PRST_attr(:),f2,'color','k','linewidth',2);
ylim([0 1])
xlim([0 1])

subplot(1,4,3)
for n = 1:length(maxtwmat)
    hold on;
    plot(PRST_attr(n,:),PRST_sim(n,:),'.','color',[1-(0.025*n),0,0.025*n],'markersize',5)
end
xlabel('Attraction')
ylabel('Similarity')
p2 = polyfit(PRST_attr(:),PRST_sim(:),1);
f2 = polyval(p2,PRST_attr(:));
hold on;plot(PRST_attr(:),f2,'color','k','linewidth',2);
ylim([0 1])
xlim([0 1])
