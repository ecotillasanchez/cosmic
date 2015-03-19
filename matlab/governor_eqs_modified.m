function [Pm_dot,P3_dot,df_dx_gov] = governor_eqs_modified(Xgov,omegas_pu,ps)
%usage:  [Pm_dot,Pi_dot,df_dx_gov] = governor_eqs_modified(Xgov,P3s,omegas_pu,ps)
% % This function outputs the right hand side of differential equations for
% the turbine-governor system model and corresponding derivatives used in jacobian.

% inputs:
%   Xgov      -> State variables of governor system
%   omegas_pu -> machines' speed in pu
%   ps        -> power system structure
%
% outputs: 
%   Pm_dot    -> rhs of differential equation
%   df_dx_gov -> df_dx for governor diff equations
%% constants
C   = psconstants;

%% governor type 1: Generic speed governing
Pm      = Xgov(1,:)';
P3      = Xgov(2,:)';
dw      = omegas_pu-1;      

% Turbine-Governor parameters
R       = ps.gov(:,C.gov.R);
Tt      = ps.gov(:,C.gov.Tt);
Ti      = ps.gov(:,C.gov.Ti);
Pref    = ps.gov(:,C.gov.Pref);
LCmax   = ps.gov(:,C.gov.LCmax);       % Upper rate (or ramp up) limit
LCmin   = ps.gov(:,C.gov.LCmin);       % Lower rate (or ramp down) limit
Pmax    = ps.gov(:,C.gov.Pmax);        % Max mechanical power limit
Pmin    = ps.gov(:,C.gov.Pmin);        % Min mechanical power limit

P3_dot            = 1./(R.*Ti) .* dw;
P1                = Pref - (dw./R + P3);
[dP2_dP1,P2]      = limiter(P1,Pmax,Pmin);        % Mechanical power limiter
Pm1_dot           = 1./Tt  .* (P2-Pm);
[dPm_dPm1,Pm_dot] = limiter(Pm1_dot,LCmax,LCmin); % rate (or ramping) limiter
    
dPm_domegas_values      = (dPm_dPm1).*(1./Tt).*(dP2_dP1).*(-1./R);
dPm_dPm_values          = (dPm_dPm1).*(-1./Tt);
dPm_dP3_values          = (dPm_dPm1).*(1./Tt).*(dP2_dP1).*(-1);
dP3_dP3_values          = zeros(length(P3),1);
dP3_domegas_values      = 1./(R.*Ti);

df_dx_gov = [dPm_domegas_values dPm_dPm_values dPm_dP3_values dP3_dP3_values dP3_domegas_values];

return;


