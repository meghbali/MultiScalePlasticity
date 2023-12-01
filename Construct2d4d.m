function [out]=Construct2d4d(coef,Imat,Hmat)

out=zeros(6);
IHmat=Hmat;
if (size(coef,2)==2)
    IHmat=Imat;
end

for i=1:size(coef,2)
    out=out+coef(i)*IHmat{i};
end

end