function [sint,cost,sinp,cosp]=tpCalc(integpoin)

cost=integpoin(4);
sint=sin(acos(cost));
if (sint==0)
    sinp=0;
    cosp=1;
else
    sinp=integpoin(3)/sint;
    cosp=integpoin(2)/sint;
end
    
end