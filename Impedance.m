function [Zlong,ZxDip,ZxQuad]=Impedance(f,gamma,sigma1,sigma2,tau1,tau2,delta,b2,b3)

mu1=1;
epsilon0=8.854187817E-12;  % Vacuum permittivity  (F/m)
mu0=pi*4E-7;    % Vacuum permeability (N/A2)
L=1;
Z0=sqrt(mu0/epsilon0);


mpole=0:1;
[~,w]=size(mpole);
% ___________________________________
b=[b2-delta b2 b3];

beta=sqrt(1-gamma^(-2));    % Beta factor
v=beta*299792458;     % Beta*Light speed   (m/s)   &&  Light speed=299792458 m/s


k=2*pi*f/v;   % Wave number

One1=f.^0;

[~,h]=size(f);
% _____________________________________________



% _____________________________________________
epsilon_c=[epsilon0*One1 ;...
    epsilon0*One1+(sigma1./(1i*2*pi*f.*(1+1i*2*pi*f*tau1))) ;...
    epsilon0*One1+(sigma2./(1i*2*pi*f.*(1+1i*2*pi*f*tau2))) ;...
    epsilon0*One1];

epsilon1=epsilon_c/epsilon0;

% _____________________________________________
nu=(abs(k)).*((1-(beta^2*epsilon1*mu1)).^0.5);

% _____________________________________________
xpp=[nu(1,:)*b(1) ; nu(2,:)*b(2) ; nu(3,:)*b(3)];    % x^{p,p}
xp1p=[nu(2,:)*b(1) ; nu(3,:)*b(2) ; nu(4,:)*b(3)];    % x^{p+1,p}
[N,~]=size(xpp);


syms XPP XP1P
digitsOld = digits(10000);
XPP=vpa(xpp);
XP1P=vpa(xp1p);
digits(digitsOld)


for r=1:w
    m=mpole(r);
    for p=1:N
        KprimeP1P=(-besselk((m-1),XP1P(p,:)))-m*(besselk(m,XP1P(p,:)))./XP1P(p,:);
        KprimePP=(-besselk((m-1),XPP(p,:)))-m*(besselk(m,XPP(p,:)))./XPP(p,:);
        IprimeP1P=(besseli((m-1),XP1P(p,:)))-m*(besseli(m,XP1P(p,:)))./XP1P(p,:);
        IprimePP=(besseli((m-1),XPP(p,:)))-m*(besseli(m,XPP(p,:)))./XPP(p,:);
        
        epsp1overnu=epsilon1(p+1,:)./nu(p+1,:);
        epspovernu=epsilon1(p,:)./nu(p,:);
        
        mup1overnu=mu1./nu(p+1,:);
        mupovernu=mu1./nu(p,:);
        
        P0=-b(p)*((nu(p+1,:)).^2./epsilon1(p+1,:));
        P111=(epsp1overnu.*besseli(m,XPP(p,:)).*KprimeP1P);
        P112=(-epspovernu.*besselk(m,XP1P(p,:)).*IprimePP);
        P121=(epsp1overnu.*besselk(m,XPP(p,:)).*KprimeP1P);
        P122=(-epspovernu.*besselk(m,XP1P(p,:)).*KprimePP);
        P211=(-epsp1overnu.*besseli(m,XPP(p,:)).*IprimeP1P);
        P212=(epspovernu.*besseli(m,XP1P(p,:)).*IprimePP);
        P221=(-epsp1overnu.*besselk(m,XPP(p,:)).*IprimeP1P);
        P222=(epspovernu.*besseli(m,XP1P(p,:)).*KprimePP);
        
        Q0=-((nu(p+1,:).^2./(nu(p,:).^2))-1).*m./(beta*epsilon1(p+1,:));
        Q11=(-besseli(m,XPP(p,:)).*besselk(m,XP1P(p,:)));
        Q12=(-besselk(m,XPP(p,:)).*besselk(m,XP1P(p,:)));
        Q21=(besseli(m,XPP(p,:)).*besseli(m,XP1P(p,:)));
        Q22=(besselk(m,XPP(p,:)).*besseli(m,XP1P(p,:)));
        
        R0=-b(p)*((nu(p+1,:)).^2./mu1);
        R111=(mup1overnu.*besseli(m,XPP(p,:)).*KprimeP1P);
        R112=(-mupovernu.*besselk(m,XP1P(p,:)).*IprimePP);
        R121=(mup1overnu.*besselk(m,XPP(p,:)).*KprimeP1P);
        R122=(-mupovernu.*besselk(m,XP1P(p,:)).*KprimePP);
        R211=(-mup1overnu.*besseli(m,XPP(p,:)).*IprimeP1P);
        R212=(mupovernu.*besseli(m,XP1P(p,:)).*IprimePP);
        R221=(-mup1overnu.*besselk(m,XPP(p,:)).*IprimeP1P);
        R222=(mupovernu.*besseli(m,XP1P(p,:)).*KprimePP);
        
        for t=1:h
            P(:,:,t)=P0(t).*[P111(t)+P112(t),P121(t)+P122(t); P211(t)+P212(t),P221(t)+P222(t)];
            Q(:,:,t)=Q0(t).*[Q11(t),Q12(t); Q21(t),Q22(t)];
            S(:,:,t)=(epsilon1(p+1,t)/mu1).*Q(:,:,t);
            R(:,:,t)=R0(t).*[R111(t)+R112(t),R121(t)+R122(t); R211(t)+R212(t),R221(t)+R222(t)];
        end
        
        M(r,p,:,:,:)=[P Q ; S R];
        
    end
    
end


% MM=zeros(w,4,4,h);
for ff=1:h
    for tt=1:w
        MM(tt,:,:,ff)=reshape(M(tt,3,:,:,ff),4,4)*reshape(M(tt,2,:,:,ff),4,4)*reshape(M(tt,1,:,:,ff),4,4);
    end
end


for vv=1:w
    for uu=1:h
        alpha(vv,uu)=(MM(vv,1,2,uu).*MM(vv,3,3,uu)-MM(vv,3,2,uu).*MM(vv,1,3,uu))./((MM(vv,1,1,uu).*MM(vv,3,3,uu)-MM(vv,1,3,uu).*MM(vv,3,1,uu)));
    end
end

Zlong=double(1i*2*pi*f*mu0*L/(2*pi*beta^2*gamma^2).*alpha(1,:));
ZxDip=double(1i*k.^2*Z0*L/(4*pi*beta*gamma^4).*alpha(2,:));
ZxQuad=double(1i*k.^2*Z0*L/(4*pi*beta*gamma^4).*alpha(1,:));
end


