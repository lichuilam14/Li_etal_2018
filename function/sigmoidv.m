

function f= sigmoidv(x,minnie,height,inflection,slope);
% f= sigmoid(x,minnie,height,inflection,slope);
% 
% x - input variable
% minnie - minimum output value
% height - height of output value
% inflection - inflection point of the sigmoide
% slope - slope of sigmoid

if nargin<2
    minne=0;
end

if nargin<1
    height=1;
end

f = minnie + height ./ (1 + exp(-(x-inflection)./slope)); %fitting function
% vals = p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));