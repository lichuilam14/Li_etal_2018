clear
clc
addpath('function');
addpath('jbtools');
Xi = 10;
fspace = -90:90;

trials = 100;
XA = randn(1,trials)*10;
XB = randn(1,trials)*10;

minval = min(XA,XB);
maxval = max(XA,XB);

for t = 1:trials;
    sigmat = 10-(normpdf(fspace,minval(t),Xi)*length(fspace));
    P = normpdf(maxval(t),fspace,sigmat);
    Rmax(t) = P*fspace';
    sigmat = 10-(normpdf(fspace,maxval(t),Xi)*length(fspace));
    P = normpdf(minval(t),fspace,sigmat);
    Rmin(t) = P*fspace';
    R(t) = Rmax(t)-Rmin(t);   
    cp(t) = sigmoidv(R(t),0,1,0,5);
    if rand<cp(t);
        ch_val(t) = maxval(t);
        uch_val(t) = minval(t);
    else
        uch_val(t) = maxval(t);
        ch_val(t) = minval(t);
    end
end


fig=figure;
set(fig,'Units','inches','position',[0,0,6,3]);
subplot(1,2,1);
plot(uch_val,-abs(R),'ko','markerfacecolor','k','markeredgecolor','w','markersize',5)
[f newx]=fitxy(uch_val,-abs(R),1);
hold on;plot(newx,f,'linewidth',2,'color','r')
xlabel('value of unchosen option')
ylabel('model output')
subplot(1,2,2);
plot(ch_val,-abs(R),'ko','markerfacecolor','k','markeredgecolor','w','markersize',5)
[f newx]=fitxy(ch_val,-abs(R),1);
hold on;plot(newx,f,'linewidth',2,'color','r')
xlabel('value of chosen option')
ylabel('model output')




