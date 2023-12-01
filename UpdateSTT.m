function [StrainMic,StressMic,StressMac,StrainMac]=UpdateSTT(ParamMicro, ...
    YieldMicro,LinkMicro,PropMicro,TensorsMicro,dx,gderiv,StrainMac, ...
    cmat,StrainMic,gparam,LoadParam,StressMic,StressMac)
%%
% updates curent micro and macro state variables


%%
% reset plastic strain to previous step
for ipoin=1:ParamMicro.npoin
    StrainMic(ipoin).plcur=StrainMic(ipoin).plpre;
end

counter1=0;
for ipoiny=1:ParamMicro.npoiny
    ipoin=YieldMicro(ipoiny).id2pt;
    for iyf=1:YieldMicro(ipoiny).id2nyf
        counter1=counter1+1;
        StrainMic(ipoin).plcur=StrainMic(ipoin).plcur+dx(counter1,1)* ...
            gderiv{ipoiny}(:,YieldMicro(ipoiny).id2yf(iyf));
    end
end

counter2=0;
for ilprog=1:6
    StrainMac.cur(ilprog,1)=StrainMac.tr(ilprog,1);
    if (LoadParam.prog(ilprog)==0)
        counter2=counter2+1;
        StrainMac.cur(ilprog,1)=StrainMac.cur(ilprog,1)+ ...
            dx(ParamMicro.nfuny+counter2,1);
    end
end

for ipoin=1:ParamMicro.npoin
    StrainMicPl{ipoin}=StrainMic(ipoin).plcur;
end

[StrainMicTemp]=StrainLocal(ParamMicro,LinkMicro,TensorsMicro, ...
    StrainMac.cur,PropMicro,StrainMicPl);

for ipoin=1:ParamMicro.npoin
    StrainMic(ipoin).cur=StrainMicTemp{ipoin};
    StressMic(ipoin).cur= ...
        cmat{ipoin}*(StrainMic(ipoin).cur-StrainMic(ipoin).plcur);
end

StrainMac.plcur=gparam.zerovec;
StressMac.cur=gparam.zerovec;
for ipoin=1:ParamMicro.npoin
    iphase=LinkMicro(ipoin).phase;
    isub=LinkMicro(ipoin).sub;
    StressMac.cur=StressMac.cur+PropMicro(iphase).frac(isub)* ...
        PropMicro(iphase).integp(isub,1)*StressMic(ipoin).cur;
    StrainMac.plcur=StrainMac.plcur+PropMicro(iphase).frac(isub)* ...
        PropMicro(iphase).integp(isub,1)*inv(ParamMicro.chom)* ...
        transpose(TensorsMicro(ipoin).A)*cmat{ipoin}*StrainMic(ipoin).plcur;
end

end