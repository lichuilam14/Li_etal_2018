%% figure 1b iii - conflict effect
clear
clc
fspace = -90:90;
Xii    = [45];  
Xjmat  = [45 -45]; 
Xsmat  = [0 5];    
maxtw  = 15;
sn     = 5;
for s = 1:length(Xsmat)
    for c = 1:2
        sigmat = maxtw -(normpdf(fspace,Xjmat(c),Xsmat(s)+sn)*length(fspace));
        sigmat(sigmat<0.1)= 0.1;
        P = normpdf(Xii,fspace,sigmat);
        R(s,c) = P*fspace';
    end
end

f = figure;
set(f,'Units','inches','position',[0,0,3,3])
bar(1,(1./R(1,1)))
hold on;bar(2,(1./R(2,1)))
hold on;bar(3,squeeze(mean(1./R(2,:),2)));
xlim([0.5 3.5])
xticks(1:3);
xticklabels({'CO','SI','RI'})
yticks([])
yticklabels({});
ylabel('simulated RT')
title('conflict - fig. 1biii')