clear 
clc
fspace=-90:90;
Xii_1 = 20;
Ximat = 1:19; % target 2
Xjmat = -45:45;
Xsmat = 20;

for j = 1:length(Xjmat)
    for i = 1:length(Ximat);
        sigmat = 10-(normpdf(fspace,Xjmat(j),Xsmat(1))*length(fspace));
        P1 = normpdf(Xii_1,fspace,sigmat);
        P2 = normpdf(Ximat(i),fspace,sigmat);
        tdiff(i) = Xii_1-Ximat(i);
        
        cp(i,j) = (sum(P2(P2>P1))+sum(P1(P1>P2)))/2;
        R1(i,j) = P1*fspace';
        R2(i,j) = P2*fspace';
        
    end
end
f = figure;
set(f,'Units','inches','position',[0,0,3,3])
imagesc(Xjmat,tdiff,cp)
set(gca,'YDir','normal')
colormap hot
ylabel('High - Low');
xlabel('Distractor orientation');
title('simple decoy - Fig.S2')