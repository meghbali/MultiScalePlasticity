function out=Convert(in,key)
%%
% conversion needed for working with vectorized form of
% tensors

if (key==1)
    factor=sqrt(2);
elseif (key==2)
    factor=1.0/sqrt(2);
end

out=in;
out(4:6,1)=out(4:6,1)*factor;

end