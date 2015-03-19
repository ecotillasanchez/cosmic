function [F, dF_dxy, dF_dxyp] = differential_algebraic_i(t,xy,xyp,ps,opt)
% usage: [F, dF_dxy, dF_dxyp] = differential_algebraic_i(t,xy,xyp,ps,opt)
% differential-algebraic implicit equations that model the elements of the power system
%
% inputs:
%   t   -> time in seconds
%   xy  -> [delta omega Pm Eap E1 Efd temp Vmag Theta]
%   xyp -> xy'
%   ps  -> power system structure
%   opt -> options structure
%
% outputs: 
%   F(1) -> ddelta_dt
%   F(2) -> domega_dt
%   F(3) -> dPm_dt
%   F(4) -> dEap_dt
%   F(5) -> dE1_dt
%   F(6) -> dEfd_dt
%   F(7) -> dtemp_dt
%   F(8) -> Pg
%   F(9) -> Qg

%% figure out some sizes
nxsmac      = 7;
n_macs      = size(ps.gen,1);
n_branches  = size(ps.branch,1);
nx          = nxsmac*n_macs;
nxy         = size(xy,1);

%% get differential and algebraic portions for F
x = xy(1:nx);
y = xy(nx+1:end);

[f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt);
[g,dg_dx,dg_dy] = algebraic_eqs(t,x,y,ps,opt);

% build output
xyp(nx+1:end)   = 0;
F               = [f;g];
F               = F - xyp; 

%% jacobian
if nargout > 1
    % build dF_dxy
    dF_dxy   = [df_dx df_dy; dg_dx dg_dy]; 
end
if nargout > 2
    % build dF_dxyp
    dF_dxyp     = sparse(1:nx,1:nx,-1,nxy,nxy);
end
