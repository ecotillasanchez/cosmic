function V = getV(ps)
% usage: V = getV(ps)
% extract the complex voltage from a power system structure

C = psconstants;
V = ps.bus(:,C.bu.Vmag) .* exp( 1i * ps.bus(:,C.bu.Vang) * pi / 180 );
