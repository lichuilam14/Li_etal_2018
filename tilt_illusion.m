%% figure 1aiii - tilt illusion simulation
clear
clc
Xi      = 0;
Xjmat   = [-45:45];
Xsmat   = [3:4:15]; 
fspace  = -90:90;
maxtw   = 10;
sn      = 6;
f = figure;
set(f,'Units','inches','position',[0,0,3,3])
for s = 1:length(Xsmat);
    for j = 1:length(Xjmat)
        sigmat = maxtw-(normpdf(fspace,Xjmat(j),Xsmat(s)+sn)*length(fspace));
        sigmat(sigmat<0.1) = 0.1;
        P      = normpdf(Xi,fspace,sigmat);
        R(s,j) = P*fspace';
    end
    hold on;
    plot(Xjmat,R(s,:),'color',[1-0.2*s,0,0+0.2*s],'linewidth',2);
end
ylim([-4 4]);
xlim([-45 45]);
line([min(Xjmat) max(Xjmat)],[0 0],'linestyle','--','color','k','linewidth',1);
ylabel('bias (deg)')
xlabel('mean flanker tilt(deg)')
title('tilt illusion - fig. 1aiii')