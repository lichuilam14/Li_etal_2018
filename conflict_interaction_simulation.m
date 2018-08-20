clear 
clc
fspace = -90:90;
Xii    = [45];   
Xjmat  = [45 -45]; 
Xsmat  = [10 20 30]; 

for s = 1:length(Xsmat)
    for c = 1:2
        sigmat = 10-(normpdf(fspace,Xjmat(c),Xsmat(s))*length(fspace));
        P = normpdf(Xii,fspace,sigmat);
        R(c,s) = P*fspace';
    end
end

f = figure;
set(f,'Units','inches','position',[0,0,3,3])
h = plot(1./R,'linewidth',2)
set(h(1),'Color',[0.6350    0.0780    0.1840]);
set(h(2),'Color',[0.4660    0.6740    0.1880]);
set(h(3),'Color',[    0    0.4470    0.7410]);
xticks([1 2])
xticklabels({'Cong','Incong'})
xlim([0.5 2.5])
title('conflict X flanker variability - fig.5c')
