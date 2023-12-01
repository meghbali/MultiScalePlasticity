function [StrainMic,StressMic,StrainMac,StressMac,ParamMicro]= ...
    Initialize(ParamMicro,gparam)


C1 = cell(ParamMicro.npoin,4);
C2 = cell(ParamMicro.npoin,2);
C3 = cell(1,6);
C4 = cell(1,4);
for k = 1:(ParamMicro.npoin*4)
    C1{k} = gparam.zerovec;
    if (k<=ParamMicro.npoin*2)
        C2{k} = gparam.zerovec;
    end
    if (k<=6)
        C3{k} = gparam.zerovec;
    end
    if (k<=4)
        C4{k} = gparam.zerovec;
    end
end

StrainMic = cell2struct(C1,{'pre','plpre','cur','plcur'},2);
StrainMac = cell2struct(C3,{'pre','plpre','cur','plcur','stage','tr'},2);
StressMic = cell2struct(C2,{'pre','cur'},2);
StressMac = cell2struct(C4,{'pre','cur','stage','tr'},2);

end