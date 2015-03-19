function [Efd_dot,E1_dot,df_dx_exc,df_dy_exc] = exciter_eqs(Xexc,mac_Vmags,ps)
% usage: [Efd_dot,E1_dot,df_dx_exc,df_dy_exc] = exciter_eqs(Xexc,mac_Vmags,ps)
% This function outputs the right hand side of differential equations for
% the exciter model and corresponding derivatives used in jacobian.

% inputs:
%   Xexc      -> State variables of exciter system
%   mac_Vmags -> machine terminal voltages in pu
%   ps        -> power system structure
%
% outputs: 
%   Efd_dot   -> rhs of differential equation
%   E1_dot    -> rhs of differential equation
%   df_dx_exc -> df_dx for exciter diff equations
%   df_dy_exc -> df_dy for exciter diff equations
%% constants
C   = psconstants;


        %% exciter type : SEXS
        Efd = Xexc(1,:)';
        E1  = Xexc(2,:)';
        
        % Exciter parameters
        Ta      = ps.exc(:,C.exc.Ta);
        Tb      = ps.exc(:,C.exc.Tb);
        Ke      = ps.exc(:,C.exc.Ke);
        Te      = ps.exc(:,C.exc.Te);
        Vref    = ps.exc(:,C.exc.Vref);
        Urmin   = ps.exc(:,C.exc.Urmin);     % Minimum value of limiter
        Urmax   = ps.exc(:,C.exc.Urmax);     % Maximum value of limiter
        
        Ve      = Vref-mac_Vmags;
        E1_dot  = 1./Tb .* (Ve-E1);
        E2      = (1-Ta./Tb).*E1+ Ve.*Ta./Tb;
        [dE3_dE2,E3]      = limiter(E2,Urmax./Ke,Urmin./Ke);    % Voltage regulator limiter
        
        % Single time constant block 
        Efd_dot = 1./Te .* (Ke.*E3-Efd);
        

        dE1_dVmag_values  = -1./Tb;
        dE1_dE1_values    = -1./Tb;
        dEfd_dVmag_values = -(Ke.*Ta.*dE3_dE2)./(Te.*Tb);
        dEfd_dE1_values = Ke.*(1-Ta./Tb).*dE3_dE2./Te;
        dEfd_dEfd_values = -1./Te;
        
        
        df_dx_exc = [dE1_dE1_values dEfd_dE1_values dEfd_dEfd_values];
        df_dy_exc = [dE1_dVmag_values dEfd_dVmag_values];


return;
