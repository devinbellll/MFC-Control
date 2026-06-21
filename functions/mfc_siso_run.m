function [yref_filter, Fnum_k, Fden_k, error, dot_dot_yref_filter, U_k, F_k] = mfc_siso_run(yref_filter_km1, yref_filter_km2, Fnum_km1, Fnum_km2, Fden_km1, Fden_km2, error_km1, error_km2, U_km1, FFilter, alpha, time, Ts, WFilter, Yref, Ym_k, pole_controller, use_control_sat)
    
% 1) Referene filter
    yref_filter = (Yref + (2*WFilter^2+2*WFilter)*yref_filter_km1+(-WFilter^2)*yref_filter_km2)/(WFilter^2+2*WFilter+1);
    %yref_filter = 90*(pi/180) + yref_filter;
    %yref_filter = Y;
    dot_yref_filter     = (yref_filter-yref_filter_km1)/Ts;
    dot_dot_yref_filter = (yref_filter-2*yref_filter_km1+yref_filter_km2)/(Ts^2);
    
% 2) Control poles (s-p)^2=s^2+as+b
    a = -2*pole_controller;
    b = -pole_controller^2;
       
    % With error
    error = Ym_k - yref_filter;
    
    % ddoty = F+alpha*u+a.dote+b.e
    sde   = -(time*error-(time-Ts)*error_km1)/Ts;
    s2d2e =  (time^2*error-2*(time-Ts)^2*error_km1+(time-2*Ts)^2*error_km2)/(Ts^2);
    sd2e  =  (time^2*error-(time-Ts)^2*error_km1)/Ts;
    de    = -time*error;
    d2e   =  time^2*error;
    d2u   =  time^2*U_km1;
    
    % d^2/ds^2 ->
    num = 2*error+4*sde+s2d2e-a*(2*de+sd2e)-b*(d2e)-alpha*d2u;  
    den = time^2;
    
% 3) Filter and estimator F    
    Fnum_k = (num+(2*FFilter^2+2*FFilter)*Fnum_km1+(-FFilter^2)*Fnum_km2)/(FFilter^2+2*FFilter+1);
    Fden_k = (den+(2*FFilter^2+2*FFilter)*Fden_km1+(-FFilter^2)*Fden_km2)/(FFilter^2+2*FFilter+1);
    
    F_k=0; 
    if (Fden_k~=0)&&(time>.1)%0.1
        F_k=Fnum_k/Fden_k;
    end

%% Command generation
    U_k = -F_k/alpha + dot_dot_yref_filter/alpha;

% Saturation
if(use_control_sat)
     U_k = max(-600, U_k);
     U_k = min(600,  U_k);
end
UFilter = 1;

%Filter of U eventually necessary
    U_k = (U_k+(UFilter-1)*U_km1)/UFilter;
 