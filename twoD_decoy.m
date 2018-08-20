%% figure 2 2D decoy effects
clear
clc
fspace = -90:90;
Xsmat  = 30;
X      = [15,10];
Y      = [10,15];
Z1     = linspace(0,30,40);
Z2     = linspace(0,30,45);
maxtw  = 10.00;
scaly  = 440;
% first attribute
for i = 1:length(Z1)
    sigmat1 = maxtw-(normpdf(fspace,Z1(i),Xsmat)*scaly);
    PX1 = normpdf(X(1),fspace,sigmat1);
    PY1 = normpdf(Y(1),fspace,sigmat1);
    X1  = PX1*fspace';
    Y1  = PY1*fspace';
    % second attribute
    for j = 1:length(Z2)
        sigmat2 = maxtw-(normpdf(fspace,Z2(j),Xsmat)*scaly);
        PX2 = normpdf(X(2),fspace,sigmat2);
        PY2 = normpdf(Y(2),fspace,sigmat2);
        X2  = PX2*fspace';
        Y2  = PY2*fspace';
        allX(i,j) = X1+X2;
        allY(i,j) = Y1+Y2;
    end
end

for n = 1:length(Z1)
    for l = 1:length(Z2)
        valdiff(n,l) = allX(n,l)-allY(n,l);
    end
end

ci = linspace(0,1,101);
for c = 1:length(ci)
colmap(c,:) = [0+ci(c),0+ci(c),1];
end
for c = 2:length(ci)
colmap(c+length(ci)-1,:) = [1,1-ci(c),1-ci(c)];
end

f = figure;
set(f,'Units','inches','position',[0,0,5,5])
contourf(Z1,Z2,valdiff',10); % z1 = x-axis 
hold on;plot(X(1),X(2),'o','markersize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','linewidth',1);
text(X(1)+1,X(2)+1,['X \in ', num2str(X)],'fontsize',12)
hold on;plot(Y(1),Y(2),'o','markersize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','linewidth',1);
text(Y(1)-2,Y(2)+2,['Y \in ',num2str(Y)],'fontsize',12)
xlabel('attribute 1 (a.u.)')
ylabel('attribute 2 (a.u.)')
colormap(colmap)
title('2D decoy effect - fig.2')
