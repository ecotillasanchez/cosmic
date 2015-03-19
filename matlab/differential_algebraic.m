function fg = differential_algebraic(t,xy,nx,ps,opt)
% usage: fg = differential_algebraic(t,xy,nx,ps,opt)

x = xy(1:nx);
y = xy((nx+1):end);

fg = [differential_eqs(t,x,y,ps,opt);
         algebraic_eqs(t,x,y,ps,opt);];
