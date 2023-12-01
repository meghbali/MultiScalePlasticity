function [YieldMicro,ParamMicro,StrainMic,ReturnIter,StrainMac,StressMic,StressMac]= ...
    ReturnMap(ParamMicro,StrainMic,LinkMicro,TensorsMicro,StrainMac, ...
    CmatLocalMicro,PropMicro,gparam,StressMac,LoadParam,StressMic)

% calculate macroscopic trial stress (prescribed strain)
StressMac.tr=ParamMicro.chom*(StrainMac.tr-StrainMac.plpre);

% get the trial micro strains
for ipoin=1:ParamMicro.npoin
    StrainMicTrPl{ipoin}=StrainMic(ipoin).plpre;
end

[StrainMicTr]=StrainLocal(ParamMicro,LinkMicro,TensorsMicro, ...
    StrainMac.tr,PropMicro,StrainMicTrPl);

%% return mapping part
ParamMicro.npoiny=0;
ParamMicro.nfuny=0;
for ipoin=1:ParamMicro.npoin
    pt2id(ipoin)=0;
end
YieldMicro(1).id2pt=1000;
YieldMicro(1).id2nyf=0;
YieldMicro(1).id2yf=[];

% calculate trial stresses and see which phases are potentially plastic
gderiv{1}=zeros(10);
for ipoin=1:ParamMicro.npoin
    iphase=LinkMicro(ipoin).phase;isub=LinkMicro(ipoin).sub;
    [~,transijkl]=RotationCalc(PropMicro(iphase).integp(isub,:),gparam);
    cmat{ipoin}=transijkl*CmatLocalMicro{iphase,isub}*transijkl';
    StressMicTr=cmat{ipoin}*(StrainMicTr{ipoin}-StrainMicTrPl{ipoin});
    
    if (strcmp(PropMicro(iphase).const,'Mohr-Coulomb')==1)
        [~,~,~,yieldreg,~]=YPEval(StressMicTr,PropMicro(iphase));
        if ((yieldreg(1)>=0) & (yieldreg(2)<=0))
            ParamMicro.npoiny=ParamMicro.npoiny+1;
            YieldMicro(ParamMicro.npoiny).id2pt=ipoin;
            pt2id(ipoin)=ParamMicro.npoiny;
            YieldMicro(ParamMicro.npoiny).id2nyf=1;
            ParamMicro.nfuny=ParamMicro.nfuny+YieldMicro(ParamMicro.npoiny).id2nyf;
            YieldMicro(ParamMicro.npoiny).id2yf=[1];
            [~,~,gderiv{ParamMicro.npoiny},~,~]=YPEval(StressMicTr,PropMicro(iphase));
        elseif ((yieldreg(1)<0) & (yieldreg(2)<0))
            ParamMicro.npoiny=ParamMicro.npoiny+1;
            YieldMicro(ParamMicro.npoiny).id2pt=ipoin;
            pt2id(ipoin)=ParamMicro.npoiny;
            YieldMicro(ParamMicro.npoiny).id2nyf=2;
            ParamMicro.nfuny=ParamMicro.nfuny+YieldMicro(ParamMicro.npoiny).id2nyf;
            YieldMicro(ParamMicro.npoiny).id2yf=[1,2];
            [~,~,gderiv{ParamMicro.npoiny},~,~]=YPEval(StressMicTr,PropMicro(iphase));
        elseif ((yieldreg(1)>0) & (yieldreg(2)>0))
            ParamMicro.npoiny=ParamMicro.npoiny+1;
            YieldMicro(ParamMicro.npoiny).id2pt=ipoin;
            pt2id(ipoin)=ParamMicro.npoiny;
            YieldMicro(ParamMicro.npoiny).id2nyf=2;
            ParamMicro.nfuny=ParamMicro.nfuny+YieldMicro(ParamMicro.npoiny).id2nyf;
            YieldMicro(ParamMicro.npoiny).id2yf=[1,3];
            [~,~,gderiv{ParamMicro.npoiny},~,~]=YPEval(StressMicTr,PropMicro(iphase));
        else
            ParamMicro.npoiny=ParamMicro.npoiny+1;
            YieldMicro(ParamMicro.npoiny).id2pt=ipoin;
            pt2id(ipoin)=ParamMicro.npoiny;
            YieldMicro(ParamMicro.npoiny).id2nyf=3;
            ParamMicro.nfuny=ParamMicro.nfuny+YieldMicro(ParamMicro.npoiny).id2nyf;
            YieldMicro(ParamMicro.npoiny).id2yf=[1,2,3];
            [~,~,gderiv{ParamMicro.npoiny},~,~]=YPEval(StressMicTr,PropMicro(iphase));
        end
    elseif (strcmp(PropMicro(iphase).const,'Drucker-Prager')==1)
        [~,yieldf,~,~,~]=YPEval(StressMicTr,PropMicro(iphase));
        if (yieldf(1)>0)
            ParamMicro.npoiny=ParamMicro.npoiny+1;
            YieldMicro(ParamMicro.npoiny).id2pt=ipoin;
            pt2id(ipoin)=ParamMicro.npoiny;
            YieldMicro(ParamMicro.npoiny).id2nyf=1;
            ParamMicro.nfuny=ParamMicro.nfuny+YieldMicro(ParamMicro.npoiny).id2nyf;
            YieldMicro(ParamMicro.npoiny).id2yf=[1];
            [~,~,gderiv{ParamMicro.npoiny},~,~]=YPEval(StressMicTr,PropMicro(iphase));
        end
    end
    
end

% update micro state variables
[StrainMic,StressMic,StressMac,StrainMac]=UpdateSTT(ParamMicro,YieldMicro,LinkMicro, ...
    PropMicro,TensorsMicro,zeros(ParamMicro.nfuny+6-sum(LoadParam.prog),1), ...
    gderiv,StrainMac,cmat,StrainMic,gparam,LoadParam,StressMic,StressMac);

% perform return mapping
kflag=0;
if (ParamMicro.npoiny==0)
    kflag=1;
end

gderiviter=0;
ReturnIter=0;
while (kflag==0)
    gderiviter=gderiviter+1;
    formatspec='----------derivative iteration= %d\n';
    fprintf(formatspec,gderiviter)
    
    %StrainMac.delcur=delstrain;
    [StressMic,StrainMic,ReturnIter,StrainMac,StressMac]= ...
        ReturnMain(TensorsMicro,YieldMicro,LinkMicro, ...
        ParamMicro,gderiv,PropMicro,StrainMic,cmat,gparam,StrainMac, ...
        LoadParam,StressMac,StressMic,pt2id);
    
    gderivpre=gderiv;
    kflag=1;
    for ipoiny=1:ParamMicro.npoiny
        ipoin=YieldMicro(ipoiny).id2pt;
        iphase=LinkMicro(ipoin).phase;isub=LinkMicro(ipoin).sub;
        
        [~,~,gderiv{ipoiny},~,~]=YPEval(StressMic(ipoin).cur,PropMicro(iphase));
        delderiv=gderiv{ipoiny}-gderivpre{ipoiny};
        for iyf=1:YieldMicro(ipoiny).id2nyf
            conflg=norm(delderiv(:,iyf))/norm(gderivpre{ipoiny}(:,iyf));
            if (conflg>1e-6)
                kflag=0;
            end
        end
    end

    % kflag=1; % avoid iteration over plastic potential derivative
    % (have negligible influence on the results)
end

% side note (important)
% phases deemed elastic at the trial step might become
% plastic after the return mapping --> the process should be re-done with
% the new number of plastic phases

end