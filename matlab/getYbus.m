function [Ybus, Yf, Yt, Yft, sh_ft] = getYbus(ps,includeShunts,useMATPOWER)
% usage: [Ybus, Yf, Yt] = getYbus(ps,includeShunts,useMATPOWER)
% Extract the Ybus matrix from power system data
% Inputs:
%  ps is the ps structure
%  inclueShunts is a bool indicating whether to put the impdeance shunts
%   into the Ybus. Default = true;
% Outputs:
%  Ybus - is the matrix
%  Yf is a matrix that calculates complex currents on the from end 
%    of each branch: If = Yf*V
%  Yt is a matrix that calculates complex currents on the to end
%    of each branch: It = Yt*V

if ~isfield(ps,'bus_i')
    ps = updateps(ps);
end

if nargin<2
    includeShunts = true;
end
if nargin>2 && useMATPOWER
    [~, bus, ~, branch] = ext2int(ps.bus, ps.gen, ps.branch);
    [Ybus,Yf,Yt] = makeYbus(ps.baseMVA, bus, branch);
    return
end
% collect some data
C = psconstants;
j = 1i;
eps = 1e-9;
n = size(ps.bus,1);
nbr = size(ps.branch,1);

% extract data from the branch matrix
F      = ps.bus_i(ps.branch(:,C.br.from));
T      = ps.bus_i(ps.branch(:,C.br.to));
status = (ps.branch(:,C.br.status)>=1);
R      = ps.branch(:,C.br.R);
X      = ps.branch(:,C.br.X);
B      = ps.branch(:,C.br.B);
shift  = ps.branch(:,C.br.shift);
tap    = ps.branch(:,C.br.tap);
tap( abs(tap)<eps ) = 1.0;

% calculate the branch impedance values
y_series = 1./(R+j*X);
tap_shift = tap .* exp(-j*pi/180 * shift);
y_tt = status.*( y_series + j*B/2);
y_ff = status.*( y_tt ./ tap.^2);
y_ft = status.*(-y_series ./ conj(tap_shift));
y_tf = status.*(-y_series ./ tap_shift);

% build these values into a sparse ybus matrix
Ybus = sparse(F,F,y_ff,n,n) + ...
       sparse(T,T,y_tt,n,n) + ...
       sparse(F,T,y_ft,n,n) + ...
       sparse(T,F,y_tf,n,n);

sh_ft = status.*(j*B/2);

% add the shunts if requested
if includeShunts
    if isfield(ps,'shunt') && ~isempty(ps.shunt)
        factor = ps.shunt(:,C.sh.frac_Y).*ps.shunt(:,C.sh.factor);
        y_shunt = factor.*(ps.shunt(:,C.sh.P) + j*ps.shunt(:,C.sh.Q))/ps.baseMVA;
        sh_bus_i = ps.bus_i(ps.shunt(:,1));
        y_sh_bus = sparse(sh_bus_i,1,y_shunt,n,1);
    else
        y_sh_bus = sparse(n,1);
    end
    Ybus = Ybus + sparse((1:n)',(1:n)',y_sh_bus,n,n);
end

% calculate Yf and Yt
if nargout > 1
    i = [(1:nbr)';(1:nbr)'];     %% double set of row indices    
    Yf = sparse(i, [F; T], [y_ff; y_ft], nbr, n);
    Yt = sparse(i, [F; T], [y_tf; y_tt], nbr, n);
    Yft = y_ft;
end

