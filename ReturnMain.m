function [StressMic,StrainMic,xiter,StrainMac,StressMac]= ...
    ReturnMain(TensorsMicro,YieldMicro,LinkMicro,ParamMicro, ...
    gderiv,PropMicro,StrainMic,cmat,gparam,StrainMac,LoadParam,StressMac, ...
    StressMic,pt2id)

%%
% implicit return mapping (n equations, n unknowns)

%%

% total number of equations and unknowns (yield function plus stress
% boundary condition
neqns=ParamMicro.nfuny+6-sum(LoadParam.prog);

% x is the solution vector, most of it is the lambda in flow rule, others
% are unknown strains due to pre-scribed stresses (so be wary of scale
% difference
dxin(1:neqns,1)=0;
dx=dxin;
dxpre=dx;

ConvFlag=0;
discard=[];
xiter=0;
disflag=0;
g=0;

% main iteration loop
while (ConvFlag==0)

    if (disflag==1)
        dx=dxin;
        dxpre=dx;
        xiter=0;
    end

    xiter=xiter+1;
    formatspec='return mapping iteration= %d\n';
    fprintf(formatspec,xiter)
    
    % initialize residual vector and tangent matrix
    ResVec=zeros(neqns,1);
    TanMat=zeros(neqns);

    % update micro state variables
    [StrainMic,StressMic,StressMac,StrainMac]=UpdateSTT(ParamMicro,YieldMicro,LinkMicro, ...
        PropMicro,TensorsMicro,dx,gderiv,StrainMac,cmat,StrainMic,gparam, ...
        LoadParam,StressMic,StressMac);
    
    % form the residual vector and tangent matrix for yielded points
    counter1=0;
    for ipoiny=1:ParamMicro.npoiny
        ipoin=YieldMicro(ipoiny).id2pt;
        iphase=LinkMicro(ipoin).phase;isub=LinkMicro(ipoin).sub;

        if (counter1==14)
            kashk=1;
        end

        [yderiv,yieldf,~,~,potf]=YPEval(StressMic(ipoin).cur,PropMicro(iphase));
        
        for iyf=1:YieldMicro(ipoiny).id2nyf
            counter1=counter1+1;

            if (ismember(counter1,discard)~=1)
                ResVec(counter1,1)=yieldf(YieldMicro(ipoiny).id2yf(iyf));
                
                counter2=0;
                for kpoiny=1:ParamMicro.npoiny
                    kpoin=YieldMicro(kpoiny).id2pt;
                    kphase=LinkMicro(kpoin).phase;ksub=LinkMicro(kpoin).sub;

                    for kyf=1:YieldMicro(kpoiny).id2nyf
                        counter2=counter2+1;

                        if (ismember(counter2,discard)~=1)
                            TanMat(counter1,counter2)= ...
                                TanMat(counter1,counter2)+ ...
                                yderiv(1:6,YieldMicro(ipoiny).id2yf(iyf))' ...
                                *PropMicro(kphase).integp(ksub,1)*cmat{ipoin}* ...
                                TensorsMicro(ipoin).D2{kpoin}* ...
                                gderiv{kpoiny}(:,YieldMicro(kpoiny).id2yf(kyf));

                            if (ipoiny==kpoiny)
                                TanMat(counter1,counter2)= ...
                                    TanMat(counter1,counter2) ...
                                    -yderiv(1:6,YieldMicro(ipoiny).id2yf(iyf))'* ...
                                    cmat{ipoin}* ...
                                    gderiv{ipoiny}(:,YieldMicro(ipoiny).id2yf(iyf));
                            end

                        end

                    end

                end
                
                interim=yderiv(1:6,YieldMicro(ipoiny).id2yf(iyf))'* ...
                        cmat{ipoin}*TensorsMicro(ipoin).A;
                counter3=0;
                for ilprog=1:6
                    if (LoadParam.prog(ilprog)==0)
                        counter3=counter3+1;
                        TanMat(counter1,counter2+counter3)= ...
                            TanMat(counter1,counter2+counter3)+interim(1,ilprog);
                    end
                end
                
            else
                TanMat(counter1,counter1)=1;
            end
            
        end
    end
    % form the residual vector and tangent matrix for yielded points
    
    % form the residual vector and tangent matrix for stress boundary
    % conditions
    counter4=0;
    for ilprog=1:6
        if (LoadParam.prog(ilprog)==0)
            counter4=counter4+1;
            ResVec(counter1+counter4,1)= ...
                StressMac.cur(ilprog,1)-StressMac.tr(ilprog,1);

            counter2=0;
            for kpoiny=1:ParamMicro.npoiny
                kpoin=YieldMicro(kpoiny).id2pt;
                kphase=LinkMicro(kpoin).phase;ksub=LinkMicro(kpoin).sub;

                for kyf=1:YieldMicro(kpoiny).id2nyf
                    counter2=counter2+1;

                    if (ismember(counter2,discard)~=1)
                        for zpoin=1:ParamMicro.npoin
                            zpoiny=pt2id(zpoin);
                            zphase=LinkMicro(zpoin).phase;zsub=LinkMicro(zpoin).sub;
                            intermilan=PropMicro(zphase).frac(zsub)*PropMicro(zphase).integp(zsub,1)* ...
                                cmat{zpoin}*PropMicro(kphase).integp(ksub,1)* ...
                                TensorsMicro(zpoin).D2{kpoin}* ...
                                gderiv{kpoiny}(:,YieldMicro(kpoiny).id2yf(kyf));
                            TanMat(counter1+counter4,counter2)=TanMat(counter1+counter4,counter2)+intermilan(ilprog,1);
                            if (zpoiny==kpoiny)
                                intermilan=PropMicro(zphase).frac(zsub)*PropMicro(zphase).integp(zsub,1)* ...
                                    cmat{zpoin}*gderiv{zpoiny}(:,YieldMicro(zpoiny).id2yf(kyf));
                                TanMat(counter1+counter4,counter2)=TanMat(counter1+counter4,counter2)-intermilan(ilprog,1);
                            end
                        end

                    end

                end

            end

            counter3=0;
            for jlprog=1:6
                if (LoadParam.prog(jlprog)==0)
                    counter3=counter3+1;
                    for kpoin=1:ParamMicro.npoin
                        kphase=LinkMicro(kpoin).phase;ksub=LinkMicro(kpoin).sub;
                        interim=cmat{kpoin}*TensorsMicro(kpoin).A;
                        TanMat(counter1+counter4,counter2+counter3)= ...
                            TanMat(counter1+counter4,counter2+counter3)+ ...
                            PropMicro(kphase).frac(ksub)* ...
                            PropMicro(kphase).integp(ksub,1)*interim(ilprog,jlprog);
                    end

                end
            end

        end
    end
    % form the residual vector and tangent matrix for stress boundary
    % conditions

    % uncommnet to check the analytical tangent matrix against the numerical one
    [TanMatNum]=PerturbTanMat(dx,StrainMac,StrainMic,StressMic,StressMac, ...
         neqns,ParamMicro,YieldMicro,LinkMicro,PropMicro,TensorsMicro,gderiv,cmat, ...
         gparam,LoadParam,discard,ResVec);

    % solve the tangential equation
    ddx=-inv(TanMat)*ResVec;
    
    dxtemp=dx+ddx;
    
    PlasConv=norm(dxtemp(1:ParamMicro.nfuny)-dxpre(1:ParamMicro.nfuny))/ ...
        norm(dxpre(1:ParamMicro.nfuny));

    counter=0;
    for ilprog=1:6
        if (LoadParam.prog(ilprog)==0)
            counter=counter+1;
            vecB(counter)=StrainMac.cur(counter,1);
            dxtempB(counter)=dxtemp(ParamMicro.nfuny+counter,1);
            dxtemppreB(counter)=dxpre(ParamMicro.nfuny+counter,1);
        end
    end
    
    if (counter==0)
        BounConv=0;
    else
        BounConv=abs((norm(dxtempB)-norm(dxtemppreB))/norm(vecB));
    end

    disflag=0;
    for ifuny=1:ParamMicro.nfuny
        if (abs(dxtemp(ifuny,1))~=dxtemp(ifuny,1))
            g=g+1;
            discard(g)=ifuny;
            disflag=1;
        end
    end

    ConvFlagY=0;
    ConvFlagB=0;
    if (disflag==0)
        
        if (ParamMicro.nfuny==0)
            ConvFlagY=1;
        elseif (norm(dxpre(1:ParamMicro.nfuny))==0)
            if (norm(dxtemp(1:ParamMicro.nfuny))==0)
                ConvFlagY=1;
            end
        elseif (PlasConv<1e-6)
            ConvFlagY=1;
        end

        if (BounConv<1e-6)
            ConvFlagB=1;
        end
    end
    
    dx=dxtemp;
    dxpre=dx;
    
    if (ConvFlagY==1 && ConvFlagB==1)
        ConvFlag=1;
    end

end
% main iteration loop

end