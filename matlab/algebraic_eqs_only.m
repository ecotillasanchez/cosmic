function [g,dg_dy] = algebraic_eqs_only(t,x,y,ps,opt)
% usage: [g,dg_dy] = algebraic_eqs_only(t,x,y,ps,opt)

%y = [y_sub;0];
if nargout>1
    [g,~,dg_dy] = algebraic_eqs(t,x,y,ps,opt);
    % drop the stuff related to delta_sys???
    %dg_dy(end,:) = [];
    %dg_dy(:,end) = [];
else
    g = algebraic_eqs(t,x,y,ps,opt);
end

% drop the stuff related to delta_sys???
%g(end) = [];
