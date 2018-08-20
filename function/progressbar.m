function out = progressbar(total,current)
clc
str    = '=';
prog   = floor((current./total)*100);
strmat = repmat(str,[1,floor(prog./5)]);
disp('running...');
disp([strmat,'>> completed ',num2str(prog),' %']);
end