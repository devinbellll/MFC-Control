function [sortie]=CSM_NG_alphaV(u1,u2,u3,u4,u5)
global paramCSM Y k fnumkm1 fnum fdenkm1 fden fconsigne fconsignekm1 erreur Ierreur;
% global Thisto Fhisto;
% global Uihisto Fihisto Yihisto; 

%%
time=u1; % pour matlab, temps
M=u2; % mesure
W=u3; % consigne
U=u4;
alpha=u5; % alpha variable
%%
%paramCSM paramètres pour la CSM
%alpha=paramCSM.alpha; 
Kp=paramCSM.Kp;
Ki=paramCSM.Ki;
UFiltre=paramCSM.UFiltre;
Umax=paramCSM.Umax;
Umin=paramCSM.Umin;
U0=paramCSM.U0;
WFiltre=paramCSM.WFiltre;
FFiltre=paramCSM.FFiltre;
y0=paramCSM.y0;
typeC=paramCSM.typeC;
Te=paramCSM.Te;

%% pahse d'initiallisation à temps=0
if time<3*Te 
    erreur=M-y0;
    Ierreur=0;
    Y=M;Ykm1=Y;dotY=0;
    U=U0;Ucsm=U0;Upid=0;k=0;
    F=0;consigne=y0;dotconsigne=0;fconsigne=y0;fconsignekm1=y0;
    fnum=0;fnumkm1=0;fden=0;fdenkm1=0;
    %disp('ini CSM')
else
%    U=U+Te;
%% Generation des signaux
% sortie et derivee
Ykm1=Y;
Y=M;
dotY=(Y-Ykm1)/Te; 
%U=u4;
      
% generation de la trajectoire de consigne et sa derivee
fconsignekm2=fconsignekm1;
fconsignekm1=fconsigne;    
fconsigne=(W+(2*WFiltre^2+2*WFiltre)*fconsignekm1+(-WFiltre^2)*fconsignekm2)/(WFiltre^2+2*WFiltre+1);
consigne=fconsigne;
if typeC==0
    fconsigne=W;
end
dotconsigne=(fconsigne-fconsignekm1)/Te;

%erreur
erreurkm1=erreur;
erreur=Y-fconsigne;
% intgrale de l'erreur pour KI (attention ajouter anti-wind-up, si
% necessaire)
Ierreur=Ierreur+(erreur+erreurkm1)/2*Te;

%% CSM nouvelle generation 
% PID
Upid=Kp*erreur+Ierreur*Ki;

%%%
% Estimation de F methode algebrique
% doty=F+alpha*u        
k=k+Te;
% -y-sdy+du
num=-Y+(k*Y-(k-Te)*Ykm1)/Te-k*alpha*U;
% 1/s^2
den=k;

% filtrage de F
fnumkm2=fnumkm1;
fnumkm1=fnum;
fdenkm2=fdenkm1;
fdenkm1=fden;
    
fnum=(num+(2*FFiltre^2+2*FFiltre)*fnumkm1+(-FFiltre^2)*fnumkm2)/(FFiltre^2+2*FFiltre+1);
fden=(den+(2*FFiltre^2+2*FFiltre)*fdenkm1+(-FFiltre^2)*fdenkm2)/(FFiltre^2+2*FFiltre+1);
%%%

F=0; % eviter la division par 0
if fden~=0
    F=fnum/fden;
end

%F=(Y-Ykm1)/Te-alpha*U;

    %F=0;
% Génération de cde
Ucsm=-F/alpha*1+Upid/alpha+dotconsigne/alpha*0;

% if time<2
%     Ucsm=U0;%disp('ici')
%     %disp(Ucsm)
% end
    

%% saturation
Ucsm=max(Umin,Ucsm);
Ucsm=min(Umax,Ucsm);

%%
%filtrage de U rarement necessaire
U=Ucsm;%(Ucsm+(UFiltre-1)*U)/UFiltre;
%disp('ici');
%disp(U)
   
end
  
sortie=[U F Upid];

