function [TanMatNum]=PerturbTanMat(dx,StrainMac,StrainMic,StressMic,StressMac, ...
    neqns,ParamMicro,YieldMicro,LinkMicro,PropMicro,TensorsMicro,gderiv,cmat, ...
    gparam,LoadParam,discard,ResVec)
%%
% calculated the tangent matrix numerically by perturbation

%%

epsil1=1e-10; % relative perturbation strength
epsil2=1e-15; % absolute perturbation strength

TanMatNum=zeros(neqns);

for indx=1:length(dx)

    StrainMacTemp=StrainMac;
    StrainMicTemp=StrainMic;
    StressMicTemp=StressMic;
    StressMacTemp=StressMac;

    dxtemp=dx;
    relative=1;
    if (dxtemp(indx,1)==0)
        dxtemp(indx,1)=dxtemp(indx,1)+epsil2;
        relative=0;
    else
        dxtemp(indx,1)=dxtemp(indx,1)*(1.0+epsil1);
    end

    [StrainMicTemp,StressMicTemp,StressMacTemp,StrainMacTemp]= ...
        UpdateSTT(ParamMicro,YieldMicro,LinkMicro,PropMicro,TensorsMicro, ...
        dxtemp,gderiv,StrainMacTemp,cmat,StrainMicTemp,gparam, ...
        LoadParam,StressMicTemp,StressMacTemp);

    ResVec2=zeros(neqns,1);

    counter1=0;
    for ipoiny=1:ParamMicro.npoiny
        ipoin=YieldMicro(ipoiny).id2pt;
        iphase=LinkMicro(ipoin).phase;isub=LinkMicro(ipoin).sub;

        [~,yieldf,~,~]=YPEval(StressMicTemp(ipoin).cur,PropMicro(iphase));

        for iyf=1:YieldMicro(ipoiny).id2nyf
            counter1=counter1+1;

            if (ismember(counter1,discard)~=1)
                ResVec2(counter1,1)=yieldf(YieldMicro(ipoiny).id2yf(iyf));
            end

        end
    end

    counter4=0;
    for ilprog=1:6
        if (LoadParam.prog(ilprog)==0)
            counter4=counter4+1;
            ResVec2(counter1+counter4,1)= ...
                StressMacTemp.cur(ilprog,1)-StressMacTemp.tr(ilprog,1);

        end
    end

    if (relative==0)
        TanMatNum(:,indx)=(ResVec2-ResVec)/epsil2;
    else
        TanMatNum(:,indx)=(ResVec2-ResVec)/(dxtemp(indx,1)*epsil1);
    end

end

for iii=1:length(discard)
    TanMatNum(iii,iii)=1;
end

end