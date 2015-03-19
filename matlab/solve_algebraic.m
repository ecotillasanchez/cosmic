function y_new = solve_algebraic(t,x,y,ps,opt)
% usage: y_new = solve_algebraic(t,x,y,ps,opt)
% solves the algebraic system (g) after a discrete change
% this assumes that ps.Ybus has been updated with the new discrete state

% set up the function to be solved
opt.nr.verbose      = opt.verbose;
opt.nr.linesearch   = opt.pf.linesearch;
g_handle = @(y) algebraic_eqs_only(t,x,y,ps,opt);
[y_new,success] = nrsolve(g_handle,y,opt);

if ~success
    y_new = [];
end
