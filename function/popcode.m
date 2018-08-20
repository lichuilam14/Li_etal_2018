function [R] =  popcode(tm,fm,fv,fspace,maxtw);
tw = (maxtw)-(normpdf(fspace,fm,fv)*length(fspace));
tw(tw<0.1)=0.1;
PR = normpdf(tm,fspace,tw);
R  = PR*fspace';