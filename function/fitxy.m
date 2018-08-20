function [f newx p s]=fitxy(x,y,polynom);
rx=ranger(x);
intv=(max(x)-min(x))./(length(x)-1);
newx=rx(1):intv:rx(2);

[p s]=polyfit(x,y,polynom);
f=polyval(p,newx);



