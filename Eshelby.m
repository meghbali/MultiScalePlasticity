function Eshelbase=Eshelby(prop,param,Imat,Hmat)

gmat=param.shear;
kmat=param.bulk;
poimat=param.pois;
emat=param.young;

for iphase=1:param.nphase
    for isub=1:prop(iphase).nsub
        aspect=prop(iphase).ar(isub);
        Eshelbase{iphase,isub}=zeros(6);
        if (strcmp(prop(iphase).type,'spherical')==1)
            cvec=[(3*kmat/(3*kmat+4*gmat)),(6*(kmat+2*gmat)/(5*(3*kmat+4*gmat)))];
        elseif (strcmp(prop(iphase).type,'long prolate')==1)
            cvec=[2/(4*(1-poimat)),2*poimat/(4*(1-poimat)),0,0,(3-4*poimat)/ ...
                (4*(1-poimat)),2*(1-poimat)/(4*(1-poimat))];
        elseif (strcmp(prop(iphase).type,'very short oblate')==1)
            p=[pi()*aspect*(3-4*poimat)-aspect^2*(16*(1-poimat)), ...
                -pi()*aspect+8*aspect^2, ...
                -pi()*aspect+8*aspect^2, ...
                8-16*poimat-2*pi()*aspect*(1-4*poimat)-32*poimat*aspect^2, ...
                pi()*(7-8*poimat)*aspect/2-16*(1-poimat)*aspect^2, ...
                8*(1-poimat)-2*pi()*aspect*(2-poimat)+8*(3-poimat)*aspect^2];
            p=p./(gmat*16*(1-poimat));
            cvec=[p(1)+2*p(2)*poimat,p(1)*poimat+p(2)*(1-poimat),p(3)+poimat*p(4), ...
                2*p(3)*poimat+(1-poimat)*p(4),p(5)*(1-2*poimat),p(6)*(1-2*poimat)];
            cvec=(2*gmat/(1-2*poimat)).*cvec;
        elseif (strcmp(prop(iphase).type,'general spheroid')==1)
            if (aspect==1)
                gamma3=1/3;
                gamma=0.5*(1-gamma3);
                phi(1)=-1/12;
                phi(2)=-1/12;
                phi(3)=phi(2);
                phi(4)=0;
                phi(5)=phi(1)/2;
                phi(6)=2*phi(2);
            else
                if (aspect>1)
                    gamma3=1/(1-aspect^2)-aspect/(1-aspect^2)* ...
                        acosh(aspect)/sqrt(aspect^2-1);
                elseif (aspect<1)
                    gamma3=1/(1-aspect^2)-aspect/(1-aspect^2)* ...
                        acos(aspect)/sqrt(1-aspect^2);
                end
                gamma=0.5*(1-gamma3);
                phi(1)=aspect^2*(4*gamma-1)-gamma;
                phi(2)=aspect^2*(1-2*gamma)-gamma;
                phi(3)=phi(2);
                phi(4)=2*(3*gamma-1);
                phi(5)=phi(1)/2;
                phi(6)=2*phi(2);
                phi=phi./(4*(1-aspect^2));
            end
            p(1)=gamma+phi(1)/(1-poimat);
            p(2)=phi(2)/(1-poimat);
            p(3)=p(2);
            p(4)=1-2*gamma+phi(4)/(1-poimat);
            p(5)=gamma+phi(1)/(2*(1-poimat));
            p(6)=(1-gamma)/2+2*phi(2)/(1-poimat);
            p=p./gmat;
            cvec=[p(1)+2*p(2)*poimat,p(1)*poimat+p(2)*(1-poimat),p(3)+poimat*p(4), ...
                2*p(3)*poimat+(1-poimat)*p(4),p(5)*(1-2*poimat),p(6)*(1-2*poimat)];
            cvec=(2*gmat/(1-2*poimat)).*cvec;
        end
        Eshelbase{iphase,isub}=Construct2d4d(cvec,Imat,Hmat);
    end
    
end
end