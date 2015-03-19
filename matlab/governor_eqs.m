function [Pm_dot,df_dx_gov] = governor_eqs(Xgov,omegas_pu,ps)
%usage:  [Pm_dot,df_dx_gov] = governor_eqs(Xgov,omegas_pu,ps)
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
dw      = omegas_pu-1;      

% Turbine-Governor parameters
R       = ps.gov(:,C.gov.R);
Tt      = ps.gov(:,C.gov.Tt);
Pref    = ps.gov(:,C.gov.Pref);
LCmax   = ps.gov(:,C.gov.LCmax);       % Upper rate (or ramp up) limit
LCmin   = ps.gov(:,C.gov.LCmin);       % Lower rate (or ramp down) limit
Pmax    = ps.gov(:,C.gov.Pmax);        % Max mechanical power limit
Pmin    = ps.gov(:,C.gov.Pmin);        % Min mechanical power limit

P1                = Pref - dw./R;
[dP2_dP1,P2]      = limiter(P1,Pmax,Pmin);        % Mechanical power limiter
Pm1_dot           = 1./Tt  .* (P2-Pm);
[dPm_dPm1,Pm_dot] = limiter(Pm1_dot,LCmax,LCmin); % rate (or ramping) limiter
    
dPm_domegas_values      = (dPm_dPm1).*(1./Tt).*(dP2_dP1).*(-1./R);
dPm_dPm_values          = (dPm_dPm1).*(-1./Tt);

df_dx_gov = [dPm_domegas_values dPm_dPm_values];

return;
