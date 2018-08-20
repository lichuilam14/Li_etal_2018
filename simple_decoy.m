%% figure 1c iii - Louie et al. distractor effect
clear
clc
fspace  = -90:90;
Xii     = [10 20]; % target vals: R1 & R2
Xjmat   = [-40:20];  % flanker means
Xjrange = max(Xjmat)-min(Xjmat);
minXj   = min(Xjmat);
Xsmat   = [10];  % flanker variances

for j = 1:length(Xjmat) % for each flanker val
    for s = 1:length(Xsmat); % for each flanker variance
        sigmat = 10-(normpdf(fspace,Xjmat(j),Xsmat(s))*length(fspace));
        P1 = normpdf(Xii(1),fspace,sigmat);
        P2 = normpdf(Xii(2),fspace,sigmat);
        P3 = normpdf(Xjmat(j),fspace,sigmat);
        cp(s,j) = (sum(P2(P2>P1))+sum(P1(P1>P2)))/2; %% what is cp?
        R1(s,j) = (P1*fspace');
        R2(s,j) = (P2*fspace');
        R3(s,j) = (P3*fspace');
    end
end

f = figure;
set(f,'Units','inches','position',[0,0,3,3])
normRdiff = (R2-R1)./R1;
normXj =  (Xjmat - minXj)./Xjrange;
normXii = (Xii - minXj)./Xjrange;
plot(normXj,normRdiff,'k','linewidth',3);
yy=get(gca,'ylim');
line([normXii(1) normXii(1)],yy,'color','b','linewidth',2);
line([normXii(2) normXii(2)],yy,'color','r','linewidth',2);
ylim([yy]);
set(gca,'box','off')
xlabel('normalised Distractor value (Z)')
ylabel('normalised difference')
title('simple decoy effect - fig.1ciii')