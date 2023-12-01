function [yideriv,yieldf,potderiv,yieldreg,potf]=YPEval(stress,propphase)
%%
% calculates the yield and plastic potential functions

%% some initial calculations
stret=VecToTen(stress);
J1=stress(1,1)+stress(2,1)+stress(3,1);

strevd=stress;
stretd=stret;
strevd(1,1)=stress(1,1)-J1/3;
strevd(2,1)=stress(2,1)-J1/3;
strevd(3,1)=stress(3,1)-J1/3;
stretd(1,1)=stret(1,1)-J1/3;
stretd(2,2)=stret(2,2)-J1/3;
stretd(3,3)=stret(3,3)-J1/3;

J2d=0.5*(strevd(1,1)^2+strevd(2,1)^2+strevd(3,1)^2+strevd(4,1)^2 ...
    +strevd(5,1)^2+strevd(6,1)^2);
sigsig=(stress(1,1)^2+stress(2,1)^2+stress(3,1)^2+stress(4,1)^2 ...
    +stress(5,1)^2+stress(6,1)^2);
J3d=det(stretd);
Jt=real(asin(-3*sqrt(3)*J3d/(2*J2d^(3/2)))/3);

c1=2*J2d^(1/2)/sqrt(3);
c2=J1/3;

%%
% principal stresses in descending order
prinstre=zeros(6,1);
prinstre(1,1)=c1*sin(Jt+2*pi()/3)+c2;
prinstre(2,1)=c1*sin(Jt)+c2;
prinstre(3,1)=c1*sin(Jt+4*pi()/3)+c2;

%% failure laws

if (strcmp(propphase.const,'Drucker-Prager')==1)

    % yield and potential, Drucker-Prager (associated)
    yieldf(1)=sqrt(3*J2d)+tan(propphase.friction)*J1/3-propphase.cohesion;
    potf(1)=yieldf(1);
    
    % yield and potential function derivatives
    potderiv=zeros(6,1);
    yideriv=zeros(6,1);
    
    c1=tan(propphase.friction)/3;
    c2=sqrt(3);

    % n1: partial J1 / partial sigma_{ij}
    n1=zeros(6,1);
    n1(1,1)=1;
    n1(2,1)=1;
    n1(3,1)=1;

    % n2: partial sqrt(J2d) / partial sigma_{ij}
    n2=strevd;
    n2=n2/(2*sqrt(J2d));
    
    yideriv=c1*n1+c2*n2;
    potderiv=yideriv;

    yieldreg(1)=0;
    yieldreg(2)=0;
    
elseif (strcmp(propphase.const,'Mohr-Coulomb')==1)
    sigy=(2*propphase.cohesion*cos(propphase.friction)) ...
        /(1-sin(propphase.friction));
    gamma=(1+sin(propphase.friction))/(1-sin(propphase.friction));

    yieldf(1)=gamma*prinstre(1,1)-prinstre(3,1)-sigy;
    yieldf(2)=gamma*prinstre(1,1)-prinstre(2,1)-sigy;
    yieldf(3)=gamma*prinstre(2,1)-prinstre(3,1)-sigy;

    potf(1)=prinstre(1,1)-prinstre(3,1);
    potf(2)=prinstre(1,1)-prinstre(2,1);
    potf(3)=prinstre(2,1)-prinstre(3,1);

    yieldreg(1)=(prinstre(1,1)-sigy/(gamma-1))/(gamma+1)-(prinstre(2,1)- ...
        sigy/(gamma-1))+(prinstre(3,1)-sigy/(gamma-1))/(gamma+1);

    yieldreg(2)=(prinstre(1,1)-sigy/(gamma-1))*gamma/(gamma+1)-(prinstre(2,1)- ...
        sigy/(gamma-1))+(prinstre(3,1)-sigy/(gamma-1))*gamma/(gamma+1);

    n1=zeros(6,1);
    n1(1,1)=1;
    n1(2,1)=1;
    n1(3,1)=1;

    n2=strevd/(2*J2d^(1/2));
    
    n3(1,1)=strevd(1,1)^2+strevd(6,1)^2/2+strevd(5,1)^2/2-2*J2d/3;
    n3(2,1)=strevd(2,1)^2+strevd(6,1)^2/2+strevd(4,1)^2/2-2*J2d/3;
    n3(3,1)=strevd(3,1)^2+strevd(4,1)^2/2+strevd(5,1)^2/2-2*J2d/3;
    n3(4,1)=strevd(5,1)*strevd(6,1)/2+strevd(2,1)*strevd(4,1)/sqrt(2)+strevd(4,1)*strevd(3,1)/sqrt(2);
    n3(5,1)=strevd(5,1)*strevd(1,1)/sqrt(2)+strevd(6,1)*strevd(4,1)/2+strevd(5,1)*strevd(3,1)/sqrt(2);
    n3(6,1)=strevd(1,1)*strevd(6,1)/sqrt(2)+strevd(2,1)*strevd(6,1)/sqrt(2)+strevd(4,1)*strevd(5,1)/2;
    n3(4,1)=n3(4,1)*sqrt(2);
    n3(5,1)=n3(5,1)*sqrt(2);
    n3(6,1)=n3(6,1)*sqrt(2);

    c(1,1)=(gamma-1)/3;
    c(1,2)=-(gamma-1)*sin(Jt)/sqrt(3)+(gamma+1)*cos(Jt)-(2/sqrt(3))*tan(3*Jt) ...
        *(-(gamma-1)*cos(Jt)/2-sqrt(3)*(gamma+1)*sin(Jt)/2);
    c(1,3)=(-1/(J2d*cos(3*Jt)))*(-(gamma-1)*cos(Jt)/2-sqrt(3)*(gamma+1)*sin(Jt)/2);

    c(2,1)=(gamma-1)/3;
    c(2,2)=-(gamma+2)*sin(Jt)/sqrt(3)+gamma*cos(Jt)-(2/sqrt(3))*tan(3*Jt) ...
        *(-(gamma/2+1)*cos(Jt)-sqrt(3)*gamma*sin(Jt)/2);
    c(2,3)=(-1/(J2d*cos(3*Jt)))*(-(gamma/2+1)*cos(Jt)-sqrt(3)*gamma*sin(Jt)/2);

    c(3,1)=(gamma-1)/3;
    c(3,2)=(2*gamma+1)*sin(Jt)/sqrt(3)+cos(Jt)-(2/sqrt(3))*tan(3*Jt) ...
        *((gamma+1/2)*cos(Jt)-sqrt(3)*sin(Jt)/2);
    c(3,3)=(-1/(J2d*cos(3*Jt)))*((gamma+1/2)*cos(Jt)-sqrt(3)*sin(Jt)/2);

    d=zeros(3);
    d(1,2)=2*cos(Jt)+2*tan(3*Jt)*sin(Jt);
    d(1,3)=(-1/(J2d*cos(3*Jt)))*(-sqrt(3)*sin(Jt));

    d(2,2)=(2/sqrt(3))*(-3*sin(Jt)/2+sqrt(3)*cos(Jt)/2)-(2/sqrt(3))*tan(3*Jt)* ...
        (-3*cos(Jt)/2-sqrt(3)*sin(Jt)/2);
    d(2,3)=(-1/(J2d*cos(3*Jt)))*(-3*cos(Jt)/2-sqrt(3)*sin(Jt)/2);

    d(3,2)=(2/sqrt(3))*(3*sin(Jt)/2+sqrt(3)*cos(Jt)/2)-(2/sqrt(3))*tan(3*Jt)* ...
        (3*cos(Jt)/2-sqrt(3)*sin(Jt)/2);
    d(3,3)=(-1/(J2d*cos(3*Jt)))*(3*cos(Jt)/2-sqrt(3)*sin(Jt)/2);

    potderiv=zeros(6,3);
    yideriv=zeros(6,3);

    for i=1:3
        yideriv(1:6,i)=c(i,1)*n1+c(i,2)*n2+c(i,3)*n3;
        potderiv(1:6,i)=d(i,1)*n1+d(i,2)*n2+d(i,3)*n3;
    end

end

end