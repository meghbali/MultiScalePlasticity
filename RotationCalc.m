function [transij,transijkl]=RotationCalc(IntegpSub,gparam)

[sint,cost,sinp,cosp]=tpCalc(IntegpSub);
transij=Rotationij(sint,cost,sinp,cosp);
transijkl=Rotationijkl(transij,gparam.guide);

end