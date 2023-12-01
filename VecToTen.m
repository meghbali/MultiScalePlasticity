function ten=VecToTen(vec)
%%
% converts a vector to tensor

vec=Convert(vec,2);

ten(1,1)=vec(1,1);
ten(2,2)=vec(2,1);
ten(3,3)=vec(3,1);
ten(2,3)=vec(4,1);
ten(1,3)=vec(5,1);
ten(1,2)=vec(6,1);
ten(2,1)=ten(1,2);
ten(3,2)=ten(2,3);
ten(3,1)=ten(1,3);

end