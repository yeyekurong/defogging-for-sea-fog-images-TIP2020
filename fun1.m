% function y=fun1(t)
% a = [-1,2,-3,9,-12,-7,-4,19];
% b = ones(1,8)*10;
% c = zeros(1,8);
% z = a + t*ones(1,8)-c;
% w = a + t*ones(1,8)-b;
% z(z>0) = 0;
% w(w<0) = 0;
% y = sum(power(z,2) + power(w,2));
% end
function y = fun1(t,L1t,channel)
a = L1t(:);
b = ones(size(a));
c = channel(:);
d = zeros(size(a));
z = a + t*b-d;
w = a + t*b-c;
z(z>0) = 0;
w(w<0) = 0;
y = sum(power(z,2) + power(w,2));
end
