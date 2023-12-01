function [Hmat]=Hillbases(gparam)

% important: when transoforming the dydic product of two tensors (which is
% a fourth-order tensor A_ijkl) into a second-order tensor, you must average two
% cases where k and l are swapped.

% identity tensor
ident=[1,0,0;0,1,0;0,0,1];

teta=[1,0,0;0,1,0;0,0,0];

% Hmat{1} to Hmat{6} are the six Hill basis tensors. These are originally
% fourth order tensors which are transformed into second order tensors for
% simplifying different operations on them
for i=1:6
    for j=1:6
        % 0.5 \teta_ij \teta_kl
        Hmat{1}(i,j)=0.5*teta(gparam.guide{i}(1,1),gparam.guide{i}(1,2))* ...
            teta(gparam.guide{j}(1,1),gparam.guide{j}(1,2));
        % \teta_ij \delta_k3 \delta_l3
        Hmat{2}(i,j)=teta(gparam.guide{i}(1,1),gparam.guide{i}(1,2))* ...
            ident(gparam.guide{j}(1,1),3)*ident(gparam.guide{j}(1,2),3);
        % \teta_kl \delta_i3 \delta_j3
        Hmat{3}(i,j)=teta(gparam.guide{j}(1,1),gparam.guide{j}(1,2))* ...
            ident(gparam.guide{i}(1,1),3)*ident(gparam.guide{i}(1,2),3);
        % \delta_i3 \delta_j3 \delta_k3 \delta_l3
        Hmat{4}(i,j)=ident(gparam.guide{i}(1,1),3)*ident(gparam.guide{i}(1,2),3)* ...
            ident(gparam.guide{j}(1,1),3)*ident(gparam.guide{j}(1,2),3);
        % 0.5*(teta_ik \teta_lj + teta_il \teta_kj - teta_ij \teta_kl)
        Hmat{5}(i,j)=0.5*(teta(gparam.guide{i}(1,1),gparam.guide{j}(1,1))*teta(gparam.guide{j}(1,2),gparam.guide{i}(1,2)) ...
            +teta(gparam.guide{i}(1,1),gparam.guide{j}(1,2))*teta(gparam.guide{j}(1,1),gparam.guide{i}(1,2)) ...
            -teta(gparam.guide{i}(1,1),gparam.guide{i}(1,2))*teta(gparam.guide{j}(1,1),gparam.guide{j}(1,2)));
        % 0.5*(\teta_ik \delta_l3 \delta_j3 + \teta_il \delta_k3 \delta_j3 + \teta_jk \delta_l3 \delta_i3 + \teta_jl \delta_k3 \delta_i3)
        Hmat{6}(i,j)=0.5*(teta(gparam.guide{i}(1,1),gparam.guide{j}(1,1))*ident(gparam.guide{j}(1,2),3)*ident(gparam.guide{i}(1,2),3) ...
            +teta(gparam.guide{i}(1,1),gparam.guide{j}(1,2))*ident(gparam.guide{j}(1,1),3)*ident(gparam.guide{i}(1,2),3) ...
            +teta(gparam.guide{i}(1,2),gparam.guide{j}(1,1))*ident(gparam.guide{j}(1,2),3)*ident(gparam.guide{i}(1,1),3) ...
            +teta(gparam.guide{i}(1,2),gparam.guide{j}(1,2))*ident(gparam.guide{j}(1,1),3)*ident(gparam.guide{i}(1,1),3));
        if (j>3)
          Hmat{1}(i,j)=Hmat{1}(i,j)*sqrt(2);
          Hmat{2}(i,j)=Hmat{2}(i,j)*sqrt(2);
          Hmat{3}(i,j)=Hmat{3}(i,j)*sqrt(2);
          Hmat{4}(i,j)=Hmat{4}(i,j)*sqrt(2);
          Hmat{5}(i,j)=Hmat{5}(i,j)*sqrt(2);
          Hmat{6}(i,j)=Hmat{6}(i,j)*sqrt(2);
        end 
    end
    if (i>3)
        Hmat{1}(i,:)=Hmat{1}(i,:)*sqrt(2);
        Hmat{2}(i,:)=Hmat{2}(i,:)*sqrt(2);
        Hmat{3}(i,:)=Hmat{3}(i,:)*sqrt(2);
        Hmat{4}(i,:)=Hmat{4}(i,:)*sqrt(2);
        Hmat{5}(i,:)=Hmat{5}(i,:)*sqrt(2);
        Hmat{6}(i,:)=Hmat{6}(i,:)*sqrt(2);
    end
end

end