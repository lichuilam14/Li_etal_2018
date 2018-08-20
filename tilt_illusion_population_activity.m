clear
clc
f = figure;
set(f,'Units','inches','position',[0,0,3,3])
fspace = linspace(-90,90,180);
Xii_1 = 15;
Xii_2 = -15;
Xii_3 = 30;
Xjmat = 45;
Xs = 10;
sigmat = 10-(normpdf(fspace,Xjmat,Xs)*length(fspace));
for n = 1:9:length(sigmat)
    tc = normpdf(fspace,fspace(n),sigmat(n));
    hold on;
    plot(fspace,tc,'color',[0.5,0.5,0.5],'linewidth',1);
end
P1 = normpdf(Xii_1,fspace,sigmat);
R1 = P1*fspace';
P2 = normpdf(Xii_2,fspace,sigmat);
R2 = P2*fspace';
P3 = normpdf(Xii_3,fspace,sigmat);
R3 = P3*fspace';
hold on;

plot(fspace,P3,'color','r','linewidth',2);
plot(fspace,P2,'color','b','linewidth',2);
plot(fspace,P1,'color',[0.5,0,0.5],'linewidth',2);

line([R3 R3],[0 0.071],'Color','r','linestyle','--','linewidth',0.8)
plot(R3,0.07,'rv','markersize',2,'linewidth',3);
plot(Xii_3,0.065,'Color',[0.5,0,0],'marker','v','markersize',2,'linewidth',3);

line([R2 R2],[0 0.045],'Color','b','linestyle','--','linewidth',0.8)
plot(R2,0.052,'bv','markersize',2,'linewidth',3);
plot(Xii_2,0.047,'Color',[0,0,0.5],'marker','v','markersize',2,'linewidth',3);

line([R1 R1],[0 0.046],'Color',[0.5,0,0.5],'linestyle','--','linewidth',0.8)
plot(R1,0.053,'color',[0.5,0,0.5],'marker','v','markersize',2,'linewidth',3);
plot(Xii_1,0.048,'Color',[0.3, 0,0.3],'marker','v','markersize',2,'linewidth',3);
xlim([-90 90])
ylabel('neural activation')
xlabel('orientation (deg)')

title('tilt illusion population activity - Fig.S1')