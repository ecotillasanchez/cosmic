function dg_dy = get_dg_dy(t,x,y,ps)

[~,~,dg_dy] = algebraic_eqs(t,x,y,ps);
