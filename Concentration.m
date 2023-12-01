function [tensors,param,CmatLocal,prop]= ...
    Concentration(param,prop,link,gparam)

% convergence tolerance for self-consistent scheme
tol=1e-06;

% construct the elasticity tensor in local coordinate for each sub-phase
for iphase=1:param.nphase
    for idir=1:prop(iphase).nsub
        if (strcmp(prop(iphase).elas,'ISO')==1)
            a1=3*prop(iphase).bulk;
            a2=2*prop(iphase).shear;
            CmatLocal{iphase,idir}=Construct2d4d([a1,a2],gparam.Imat,gparam.Hmat);
        elseif (strcmp(prop(iphase).elas,'TISO')==1)
            dd=prop(iphase).young^2*prop(iphase).young3^2/((1+prop(iphase).pois)* ...
                ((1-prop(iphase).pois)*prop(iphase).young3-2*prop(iphase).pois3^2* ...
                prop(iphase).young));
            c11=dd*(1/(prop(iphase).young*prop(iphase).young3)-prop(iphase).pois3^2/ ...
                prop(iphase).young3^2);
            c12=dd*(prop(iphase).pois/(prop(iphase).young*prop(iphase).young3)+prop(iphase).pois3^2/ ...
                prop(iphase).young3^2);
            c13=dd*(1+prop(iphase).pois)*prop(iphase).pois3/(prop(iphase).young* ...
                prop(iphase).young3);
            c33=dd*(1-prop(iphase).pois^2)/prop(iphase).young^2;
            b5=prop(iphase).shear3;
            b4=(c11-c12)/2;
            b3=c33;
            b2=c13;
            b1=(c11+c12)/2;
            CmatLocal{iphase,idir}=Construct2d4d([2*b1,b2,b2,b3,2*b4,2*b5],gparam.Imat,gparam.Hmat);
        elseif (strcmp(prop(iphase).elas,'custom')==1)
            CmatLocal{iphase,idir}=Construct2d4d([prop(iphase).vec(1,idir), ...
                prop(iphase).vec(2,idir),prop(iphase).vec(3,idir), ...
                prop(iphase).vec(4,idir),prop(iphase).vec(5,idir), ...
                prop(iphase).vec(6,idir)],gparam.Imat,gparam.Hmat);
        end
    end
end

%%% Elastic tensor and concentrations tensors
param.bulk=param.inbulk;
param.young=param.inyoung;
param.pois=param.inpois;
param.shear=param.inshear;

a1=3*param.bulk;
a2=2*param.shear;
param.chom=Construct2d4d([a1,a2],gparam.Imat,gparam.Hmat);

selflag=1;
iter=0;
while (selflag==1)
    iter=iter+1;
    
    % calculate the base Eshelby tensor in local coordinate
    Eshelbase=Eshelby(prop,param,gparam.Imat,gparam.Hmat);
    
    % FA == \sum \phi A
    FA=gparam.zeroten;
    % FAP == \sum \phi A P
    FAP=gparam.zeroten;
    for ipoin=1:param.npoin
        iphase=link(ipoin).phase;
        isub=link(ipoin).sub;
        [sint,cost,sinp,cosp]=tpCalc(prop(iphase).integp(isub,:));
        
        pangle=acos(cosp);
        if (sinp<0)
            pangle=2*pi()-pangle;
        end
        prop(iphase).angle(isub,1:2)=[acos(cost),pangle];
        transij=Rotationij(sint,cost,sinp,cosp);
        transijkl=Rotationijkl(transij,gparam.guide);
        
        % rotate local to global
        Eshel=transijkl*Eshelbase{iphase,isub}*transijkl';
        cmat{ipoin}=transijkl*CmatLocal{iphase,isub}*transijkl';
        Ptensor{ipoin}=Eshel*inv(param.chom);
        Hillinv=gparam.I0mat+Ptensor{ipoin}*(cmat{ipoin}-param.chom);
        Hill{ipoin}=inv(Hillinv);
        % directional integration OR summation
        FA=FA+prop(iphase).frac(isub)*Hill{ipoin}* ...
            prop(iphase).integp(isub,1);
        FAP=FAP+prop(iphase).frac(isub)*Hill{ipoin}*Ptensor{ipoin}* ...
            prop(iphase).integp(isub,1);
    end
    
    param.chom=gparam.zeroten;
    for ipoin=1:param.npoin
        iphase=link(ipoin).phase;
        isub=link(ipoin).sub;
        
        tensors(ipoin).A=Hill{ipoin}*inv(FA);
        param.chom=param.chom+prop(iphase).frac(isub)*cmat{ipoin}* ...
            tensors(ipoin).A*prop(iphase).integp(isub,1);
    end
    
    ghom=param.chom(4,4)/2;
    khom=(param.chom(1,2)*3+2*ghom)/3;
    younghom=9*khom*ghom/(3*khom+ghom);
    poihom=(3*khom-2*ghom)/(6*khom+2*ghom);
    
    if (strcmp(param.scheme,'mori-tanaka')==1)
        selflag=0;
        param.bulk=khom;
        param.shear=ghom;
        param.pois=poihom;
        param.young=younghom;
    elseif (strcmp(param.scheme,'self-consistent')==1)
        if ((abs((khom-param.bulk)/param.bulk)<tol) & (abs((ghom-param.shear)/param.shear)<tol))
            selflag=0;
            formatspec='self-consistent iterations= %d\n';
            fprintf(formatspec,iter)
        end
        param.bulk=khom;
        param.shear=ghom;
        param.pois=poihom;
        param.young=younghom;
    end

end

%%% Influence tensors
FCAP=zeros(6);
for ipoin=1:param.npoin
    iphase=link(ipoin).phase;
    isub=link(ipoin).sub;
    FCAP=FCAP+prop(iphase).frac(isub)*(param.chom-cmat{ipoin})*Hill{ipoin}* ...
        Ptensor{ipoin}*prop(iphase).integp(isub,1);
end

for i=1:6
    if (FCAP(i,i)==0)
        FCAP(i,i)=1e-15;
    end
end

for ipoin=1:param.npoin

    iphase=link(ipoin).phase;
    isub=link(ipoin).sub;
    itensor=tensors(ipoin).A;
    for jpoin=1:param.npoin
        jphase=link(jpoin).phase;
        jsub=link(jpoin).sub;
        jtensor=tensors(jpoin).A;
        HH=-itensor*Hill{jpoin}*Ptensor{jpoin}+(itensor*FAP- ...
            Hill{ipoin}*Ptensor{ipoin})*inv(FCAP)*((param.chom-cmat{jpoin}) ...
            *Hill{jpoin}*Ptensor{jpoin}+transpose(gparam.I0mat-jtensor));
        
        % D is used when working with pre-stress
        tensors(ipoin).D{jpoin}=prop(jphase).frac(jsub)*HH;
        % D2 is used when working with pre-strain
        tensors(ipoin).D2{jpoin}=prop(jphase).frac(jsub)*HH*cmat{jpoin};
        
        if (ipoin==jpoin)
            tensors(ipoin).D{jpoin}=gparam.I0mat*Hill{ipoin}*Ptensor{ipoin}/prop(jphase).integp(jsub,1)+ ...
                prop(jphase).frac(jsub)*HH;
            
            tensors(ipoin).D2{jpoin}=(gparam.I0mat*Hill{ipoin}*Ptensor{ipoin}/prop(jphase).integp(jsub,1))*cmat{jpoin}+ ...
                prop(jphase).frac(jsub)*HH*cmat{jpoin};
        end
    end
end