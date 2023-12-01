function [trans3d]=Rotationij(sint,cost,sinp,cosp)
%%
% rotatation matrix for a second order tensor

%%
% rotation p about x3
trans1=[cosp,sinp,0;-sinp,cosp,0;0,0,1];

% rotation t about x2
trans2=[cost,0,-sint;0,1,0;sint,0,cost];

% % rotation t about x1
% trans2=[1,0,0;0,cost,-sint;0,sint,cost];

trans3d=trans2*trans1;
trans3d=inv(trans3d);
end