function [StrainMac]=SolveMacroStrain(ParamMicro,StrainMac,StressMac,LoadParam)
%% 
% Calculates full trial strain tensor

%%
chom=ParamMicro.chom;
strainel=StrainMac.cur-StrainMac.plpre;
stress=StressMac.cur;
for idx=1:6
    if (LoadParam.prog(idx)==1)
        stress=stress-ParamMicro.chom(:,idx)*strainel(idx,1);
        chom(idx,:)=0;chom(:,idx)=0;chom(idx,idx)=1;
    end
end
for idx=1:6
    if (LoadParam.prog(idx)==1)
        stress(idx,1)=strainel(idx,1);
    end
end
StrainMac.tr=chom\stress+StrainMac.plpre;
StrainMac.cur=StrainMac.tr;

end