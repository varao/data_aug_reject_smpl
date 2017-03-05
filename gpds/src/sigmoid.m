function y = sigmoid(x,offset)
   y = 1./(1+exp(-1.*(x-offset)));

% y = x;
% y(y>1) = 1;
% y(y<0) = 0.01;
