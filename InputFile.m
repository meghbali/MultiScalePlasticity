function [LinkMicro,LoadParam,gparam,ParamMicro,PropMicro]=InputFile()
%% 
% define the material properties and set/load some general parameters

%% DEFINE SOME PARAMETERS & AUXILARY VECTORS & TENSORS

% the following shows how the 3X3 representation (tensor) is mapped to a 
% 1X6 (vector) representation. Basically, how the component ij of a 
% second order tensor is mapped to component i of the corresponding vector 
gparam.guide{1}=[1,1];
gparam.guide{2}=[2,2];
gparam.guide{3}=[3,3];
gparam.guide{4}=[2,3];
gparam.guide{5}=[1,3];
gparam.guide{6}=[1,2];

% load Hill basis functions
[gparam.Hmat]=Hillbases(gparam);

%Imat: 0.5*(\delta_ik \delta_jl + \del_il \del_jk)
gparam.I0mat=gparam.Hmat{1}+gparam.Hmat{4}+gparam.Hmat{5}+gparam.Hmat{6};
%Jmat: (del_ij \del_kl)/3
gparam.Imat{1}=(2*gparam.Hmat{1}+gparam.Hmat{2}+gparam.Hmat{3}+ ...
    gparam.Hmat{4})/3;
%Kmat
gparam.Imat{2}=(gparam.Hmat{1}-gparam.Hmat{2}-gparam.Hmat{3}+ ...
    2*gparam.Hmat{4}+3*gparam.Hmat{5}+3*gparam.Hmat{6})/3;
%Kmat check
% Imat{2}=I0mat-Imat{1};

gparam.nnvec=zeros(6,1);gparam.nnvec(3,1)=1;
gparam.eyevec=zeros(6,1);gparam.eyevec(1,1)=1;gparam.eyevec(2,1)=1;gparam.eyevec(3,1)=1;
gparam.zerovec=zeros(6,1);
gparam.zeroten=zeros(6);

% numerical integration points
r1=0.5;
s1=(sqrt(5)+1)/4;
t1=(sqrt(5)-1)/4;
a1=1/15;
r2=sqrt((9-4*sqrt(3))/33);
s2=sqrt((15+8*sqrt(3))/33);
t2=sqrt(1/3);
u2=sqrt((15-8*sqrt(3))/33);
v2=sqrt((9+4*sqrt(3))/33);
a2=9/280;
b2=(122+9*sqrt(3))/3360;
c2=(122-9*sqrt(3))/3360;
stroud1=[a1,r1,s1,t1;a1,r1,-s1,t1;a1,-r1,s1,t1;a1,-r1,-s1,t1; ...
    a1,t1,r1,s1;a1,t1,-r1,s1;a1,-t1,r1,s1;a1,-t1,-r1,s1;a1,s1,t1,r1; ...
    a1,s1,-t1,r1;a1,-s1,t1,r1;a1,-s1,-t1,r1;a1,1,0,0;a1,0,1,0;a1,0,0,1];
stroud2=[a2,t2,t2,t2;a2,t2,-t2,t2;a2,-t2,t2,t2;a2,-t2,-t2,t2; ...
    b2,s2,r2,r2;b2,s2,-r2,r2;b2,-s2,r2,r2;b2,-s2,-r2,r2;b2,r2,s2,r2; ...
    b2,r2,-s2,r2;b2,-r2,s2,r2;b2,-r2,-s2,r2;b2,r2,r2,s2;b2,r2,-r2,s2; ...
    c2,-u2,-v2,v2;c2,v2,u2,v2;c2,v2,-u2,v2;c2,-v2,u2,v2;c2,-v2,-u2,v2; ...
    c2,v2,v2,u2;c2,v2,-v2,u2;c2,-v2,v2,u2;c2,-v2,-v2,u2;b2,-r2,r2,s2; ...
    b2,-r2,-r2,s2;c2,u2,v2,v2;c2,u2,-v2,v2;c2,-u2,v2,v2];

%% MICRO PARAMETERS

%%% geometrical features

% rev volume
revsize=(4/3)*pi()*1^3;

% number of phases
ParamMicro.nphase=2;

% phase shapes
% the 'general spheroid' can be used instead of 'long prolate' and
% 'very short oblate' given that appropriate aspect ratios are assigned
% but it doesn't work for aspect ratio=1. you should use aspect
% ratio=0.9999.
% 'general spheroid', 'spherical', 'long prolate', 'very short oblate'
PropMicro(1).type='general spheroid';
PropMicro(1).ar(1)=1000;
PropMicro(2).type='spherical';
PropMicro(2).ar(1)=1;

% directional discretization
PropMicro(1).tetn=1;
PropMicro(2).tetn=1;
PropMicro(1).tetp=1;
PropMicro(2).tetp=1;
% total number of sub-phases
ninc1=(PropMicro(1).tetn-1)*PropMicro(1).tetp+1;
ninc2=(PropMicro(2).tetn-1)*PropMicro(2).tetp+1;

% number of phases in each teta family
PropMicro(1).nuteta=1;
PropMicro(2).nuteta=1;
% number density of phases in each teta family (1/mm^3)
PropMicro(1).nubarteta=PropMicro(2).nuteta/revsize;
PropMicro(2).nubarteta=PropMicro(2).nuteta/revsize;
% initial number density of phases in each teta family (this can optionally change 
% as the REV changes volume during loading)
PropMicro(1).nubartetain=PropMicro(2).nubarteta;
PropMicro(2).nubartetain=PropMicro(2).nubarteta;

% Integration angles for discrete sub-phases
kounter=0;
for ifam=1:(PropMicro(2).tetn+1)
    if (ifam==1 | ifam==(PropMicro(2).tetn+1))
        jmax=1;
    else
        jmax=PropMicro(2).tetp;
    end
    for jfam=1:jmax
        kounter=kounter+1;
        uniform=(ifam-1)/PropMicro(2).tetn;
        t(kounter)=acos(1-2*uniform);
        t(kounter)=(ifam-1)*pi()/PropMicro(2).tetn;
        p(kounter)=(jfam-1)*2*pi()/PropMicro(2).tetp;
        discret(kounter,4)=cos(t(kounter));
        discret(kounter,3)=sin(t(kounter))*sin(p(kounter));
        discret(kounter,2)=sin(t(kounter))*cos(p(kounter));
        discret(kounter,1)=1;
    end
end

% phase volume fraction
PropMicro(2).frac(1)=0.32; %PropMicro(2).nubarteta*(4/3)*pi()*radi2^3*asra2
PropMicro(1).frac(1)=0.68; %1-ninc2*PropMicro(2).frac(1);
% initial phase volume fraction (can optionally change as the REV changes
% volume during loading)
PropMicro(1).fracin(1)=PropMicro(1).frac(1);
PropMicro(2).fracin(1)=PropMicro(2).frac(1);

% which integration scheme to be used: stroud1, stroud2 OR custom
% [1,5,0,1] or [1,1,0,0]
PropMicro(1).integp=stroud1;
PropMicro(2).integp=[1,1,0,1];

% for isub=1:size(discret,1)
%     PropMicro(2).integp(isub,:)=discret(isub,:);
% end

% find number of sub-phases
PropMicro(1).nsub=size(PropMicro(1).integp,1);
PropMicro(2).nsub=size(PropMicro(2).integp,1);

ParamMicro.npoin=0;
for iphase=1:ParamMicro.nphase
    for isub=1:PropMicro(iphase).nsub
        ParamMicro.npoin=ParamMicro.npoin+1;
        PropMicro(iphase).ar(isub)=PropMicro(iphase).ar(1);
        PropMicro(iphase).arin(isub)=PropMicro(iphase).ar(isub);
        PropMicro(iphase).frac(isub)=PropMicro(iphase).frac(1);
        PropMicro(iphase).fracin(isub)=PropMicro(iphase).fracin(1);
        LinkMicro(ParamMicro.npoin).phase=iphase;
        LinkMicro(ParamMicro.npoin).sub=isub;
    end
end

%%% mechanical properties of phases

% constitutive law: Elastic, Drucker-Prager, Mohr-Coulomb
PropMicro(1).const='Drucker-Prager';
PropMicro(2).const='Elastic';

% elasticity model: isotropic: ISO, transversely isotropic: TISO, 
% custom elasticity: custom
PropMicro(1).elas='ISO';
PropMicro(2).elas='ISO';

% Young's modulus
PropMicro(1).young=114037; % MPa
PropMicro(2).young=0;
% Poisson'a ratio
PropMicro(1).pois=0.2699;
PropMicro(2).pois=0.5;
% Bulk modulus
PropMicro(1).bulk=PropMicro(1).young/(3*(1-2*PropMicro(1).pois));
PropMicro(2).bulk=2300; %PropMicro(2).young/(3*(1-2*PropMicro(2).pois));
% Shear modulus
PropMicro(1).shear=PropMicro(1).young/(2*(1+PropMicro(1).pois));
PropMicro(2).shear=PropMicro(2).young/(2*(1+PropMicro(2).pois));

% define plasticity parameters
% friction angle (radian)
PropMicro(1).friction=57.8*pi()/180;
% cohesion (MPa)
PropMicro(1).cohesion=82.3;

%% HOMOGENIZATION SCHEME
% homogenization scheme to be used: currently only 'mori-tanaka' and 
% 'self-consistent' are available
ParamMicro.scheme='self-consistent';

% the phase to be considered matrix
ParamMicro.imat=1;
% initial guesses (self-consistent scheme)
ParamMicro.inbulk=PropMicro(ParamMicro.imat).bulk;
ParamMicro.inshear=PropMicro(ParamMicro.imat).shear;
ParamMicro.inyoung=PropMicro(ParamMicro.imat).young;
ParamMicro.inpois=PropMicro(ParamMicro.imat).pois;

%% LOADING PROGRAM
% number of loading steps in each loading stage
LoadParam(1).nstep=0;
LoadParam(2).nstep=100;
% type of loading: 0 means stress-controled, 1 means strain controlled
LoadParam(1).prog=[0,0,0,0,0,0]';
LoadParam(2).prog=[0,0,1,0,0,0]';
% strain loading
LoadParam(1).StrainStep=[0,0,0,0,0,0]'/LoadParam(1).nstep;
LoadParam(2).StrainStep=[0,0,0.002,0,0,0]'/LoadParam(2).nstep;
% stress loading
LoadParam(1).StressStep=[0,0,0,0,0,0]'/LoadParam(1).nstep;
LoadParam(2).StressStep=[0,0,0,0,0,0]'/LoadParam(2).nstep;

%% NOTE
% for a phase with sub-phases: if the phase is discrete and you want each 
% sub-phase be counted as a separate phase, then assign the sub-phases 
% fractions to the whole phase and make the weights for each sub-phase one.
% if the phase is continuous and you want the sub-phases be like integration
% points, assign the total fraction to the phase and use the stroud's
% integration weights.
% in the second case you shouldn't treat the integration points as phases.
% for instance, the strain induced in another phase due to plastic strain
% in one of the integration points doesn't have any physical meaning. 
% It turns out the influence tensors in two cases are different in a 
% multiplication by weighting factor (continuous is bigger). 
% However, the opposite case still makes physical sense and is equal to the 
% corresponsing value in the dicrete case.
end