function [strainlocal]=StrainLocal(param,link,tensors,strain,prop,strainpl)
% caculates the localized strains

for ipoin=1:param.npoin
    iphase=link(ipoin).phase;isub=link(ipoin).sub;
    itensor=tensors(ipoin).A;
    strainlocal{ipoin}=itensor*strain;
    for jpoin=1:param.npoin
        jphase=link(jpoin).phase;jsub=link(jpoin).sub;
        strainlocal{ipoin}=strainlocal{ipoin}+prop(jphase).integp(jsub,1)* ...
            tensors(ipoin).D2{jpoin}*strainpl{jpoin};
    end
end