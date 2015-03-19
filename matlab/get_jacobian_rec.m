function [J,df_dx,df_dy,dg_dx,dg_dy] = get_jacobian_rec(t,xy,nx,ps,opt)
% usage: [J,df_dx,df_dy,dg_dx,dg_dy] = get_jacobian_rec(t,xy,nx,ps,opt)

x = xy(1:nx);
y = xy((1+nx):end);

[~, df_dx,df_dy] = differential_eqs_rec(t,x,y,ps,opt);
[~, dg_dx,dg_dy] = algebraic_eqs_rec(t,x,y,ps,opt);

J = [df_dx df_dy; dg_dx dg_dy];
