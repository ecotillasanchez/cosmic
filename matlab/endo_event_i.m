function [value,isterminal,direction] = endo_event_i(~,xy,~,ix,ps)
% usage: [value,isterminal,direction] = endo_event_i(~,xy,~,ix,ps)
% defines an exogenous event (relay) that will stop the DAE integration

% constants and settings
C = psconstants;
oc_threshold            = ps.relay(ix.re.oc,C.re.threshold);   % over-current threshold
Vmag_threshold          = ps.relay(ix.re.uvls,C.re.threshold); % undervoltage threshold
omega_pu_threshold      = ps.relay(ix.re.ufls,C.re.threshold); % underfrequency threshold
dist_zone1_threshold    = ps.relay(ix.re.dist,C.re.threshold); % distance relay zone 1 setting

% extract some info from the inputs
x               = xy(1:ix.nx); 
y               = xy(ix.nx+1:end);
Vmags           = y(ix.y.Vmag);
Thetas          = y(ix.y.theta);
V               = Vmags .* exp(1i*Thetas);
If              = ps.Yf * V;
Imag_f          = abs(If);

% temp_threshold  = ps.relay(ix.re.temp,C.re.threshold);
F               = ps.bus_i(ps.branch(:,C.br.from));
y_apparent      = Imag_f./Vmags(F);
nload           = size(ps.shunt,1);
load_freq       = zeros(nload,1);
near_gen        = ps.shunt(:,C.sh.near_gen);
for i = 1:nload
    near_gen_id = ps.gen_i(near_gen(i));
    if near_gen_id ~= 0
        load_freq_source = ix.x.omega_pu(near_gen_id);
        load_freq(i)  = x(load_freq_source);
    end
end

% get the shunt_bus locations
sh_bus_ix = ps.bus_i(ps.shunt(:,1));
Vmag_sh   = y(ix.y.Vmag(sh_bus_ix));
% build zero crossing function
value = [oc_threshold - Imag_f;                     % trigger over current relay
         Vmag_sh              - Vmag_threshold;     % trigger undervoltage load shedding
         load_freq            - omega_pu_threshold; % trigger underfrequency load shedding
         dist_zone1_threshold - y_apparent];        % trigger zone 1 distance relay
     
% for relays that have already tripped set value to be in a safe range
is_tripped = ps.relay(:,C.re.tripped)==1;
value(is_tripped) = 10;
%value(ix.re.ufls) = 10;
%value(ix.re.dist) = 10;

% set directions
n_relays = size(ps.relay,1);
isterminal  =    ones(n_relays,1);    % terminate integration
direction   =   -ones(n_relays,1);    % to detect crossing from + to -
