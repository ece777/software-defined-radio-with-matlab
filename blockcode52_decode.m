% blockcode52.m Part 1: Definition of (5,2) binary linear block code% the generator and parity check matrices
function z=blockcode52_decode(y);
% blockcode52.m Part 2: encoding and decoding data
%m=10000;                         % length of message
%dat=0.5*(sign(rand(1,m)-0.5)+1); % m random 0s and 1s
m = length(y);
mod_dat = mod(m, 10);
nmod = 0;
if(mod_dat ~= 0)
    nmod = 10 - mod_dat;
end
 
for(j=1:nmod)
    y(m+j) = 0;
end
 
p=.0;                            % probability of bit flip
m = length(y)/2.5;
 
%Decoding
k1=1;m
for i=1:2:m
  y1(1) = y(k1);
  y1(2) = y(k1+1);
  y1(3) = y(k1+2);
  y1(4) = y(k1+3);
  y1(5) = y(k1+4);  
  k1 = k1+5;
  eh=mod(y1*h',2);                % multiply by parity check h'
  ehind=eh(1)*4+eh(2)*2+eh(3)+1;  % turn syndrome into index
  e=syn(ehind,:);                 % error from syndrome table
  y1=mod(y1-e,2);                 % add e to correct errors
  for j=1:max(size(x))            % recover message from codewords
    if y1==cw(j,:), z(i:i+1)=x(j,:); end
  end
end
