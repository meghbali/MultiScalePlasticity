function transijkl=Rotationijkl(transij,guide)
%%
% rotatation matrix for a fourth order tensor

%%
for i=1:6
    for j=1:6
        transten1(i,j)=transij(guide{i}(1,1),guide{j}(1,1))* ...
            transij(guide{i}(1,2),guide{j}(1,2));
        transten2(i,j)=transij(guide{i}(1,1),guide{j}(1,2))* ...
            transij(guide{i}(1,2),guide{j}(1,1));
        transijkl(i,j)=0.5*(transten1(i,j)+transten2(i,j));
        if (j>3)
          transijkl(i,j)=transijkl(i,j)*sqrt(2);
        end
    end
    if (i>3)
        transijkl(i,:)=transijkl(i,:)*sqrt(2);
    end
end

end