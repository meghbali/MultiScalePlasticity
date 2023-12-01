clear
clc

load("/MATLAB Drive/AllVars.mat")

for i=1:(100+1)
    strain1(i)=out(i).strain(3);
    stress1(i)=out(i).stress(3);
end

plot(strain1,stress1,'--r')

% hold on
% 
% load("/MATLAB Drive/New.mat")
% 
% for i=1:(100+1)
%     strain2(i)=out(i).strain(3);
%     stress2(i)=out(i).stress(3);
% end
% 
% plot(strain2,stress2,'-.b')