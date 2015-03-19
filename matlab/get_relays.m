function relay = get_relays(ps,mode,opt)
% usage: relay = get_relays(ps,mode,opt)
% gets the relay matrix for the system described in ps

C = psconstants;
m = size(ps.branch,1);
n_shunt = size(ps.shunt,1);
% n_macs = size(ps.gen,1);
if nargin<3, opt=psoptions; end
if nargin<2, mode='all'; end


switch mode
    case 'temperature'
        relay = zeros(m,C.relay.cols);
        br_status = ps.branch(:,C.br.status);
        relay(:,C.re.type)              = C.relay.temp;
        relay(:,C.re.setting1)          = ps.branch(:,C.br.rateB)./ps.branch(:,C.br.rateA);
        relay(:,C.re.setting2)          = ps.branch(:,C.br.rateC)./ps.branch(:,C.br.rateA);
        relay(:,C.re.threshold)         = relay(:,C.re.setting1).^2.*(opt.sim.temp.TA);
        relay(:,C.re.tripped)           = ~br_status;
        relay(:,C.re.branch_loc)        = ps.branch(:,C.br.id);
        relay(:,C.re.temp_K)            = log((relay(:,C.re.threshold)./opt.sim.temp.TA-relay(:,C.re.setting2).^2)./(1-relay(:,C.re.setting2).^2))/(-60); % T_ambient = 20 degree C
        relay(:,C.re.temp_R)            = opt.sim.temp.TA./(ps.branch(:,C.br.rateA)/ps.baseMVA).^2.*relay(:,C.re.temp_K); 
        relay(:,C.re.state_a)           = relay(:,C.re.temp_R).*ps.branch(:,C.br.Imag_f).^2./relay(:,C.re.temp_K);
        relay(:,C.re.id)                = 1:m;
    case 'oc'   % over current relay
        relay = zeros(m,C.relay.cols);
        br_status = ps.branch(:,C.br.status);
        relay(:,C.re.type)              = C.relay.oc;
        Imax            = ps.branch(:,C.br.rateB)/ps.baseMVA;
        overload_max    = Imax * 1.5 * 5;    % 5 seconds at a 50% overload
        relay(:,C.re.setting1)          = Imax;
        relay(:,C.re.threshold)         = overload_max;
        relay(:,C.re.tripped)           = ~br_status;
        relay(:,C.re.branch_loc)        = ps.branch(:,C.br.id);
        relay(:,C.re.id)                = 1:m;
    case 'uvls' % under voltage load shedding relay settings
        relay = zeros(n_shunt,C.relay.cols);
        relay(:,C.re.type)      = C.relay.uvls;
        relay(:,C.re.shunt_loc) = ps.shunt(:,C.sh.bus);
        relay(:,C.re.threshold) = opt.sim.uvls_limit;
        relay(:,C.re.setting1)  = opt.sim.uvls_delta; % the amount of load to shedd on undervoltage
        relay(:,C.re.tripped)   = 0;
        relay(:,C.re.id)        = 1:n_shunt;
    case 'ufls'
        relay = zeros(n_shunt,C.relay.cols);
        relay(:,C.re.type) = C.relay.ufls;
        relay(:,C.re.shunt_loc) = ps.shunt(:,C.sh.bus);
        relay(:,C.re.threshold) = opt.sim.ufls_limit;
        relay(:,C.re.setting1) = opt.sim.ufls_delta;
        relay(:,C.re.tripped)   = 0;
        relay(:,C.re.id)        = 1:n_shunt;

    case 'distance'
        relay = zeros(m,C.relay.cols);
        relay(:,C.re.type) = C.relay.dist;
        relay(:,C.re.branch_loc) = ps.branch(:,C.br.id);
        R = ps.branch(:,C.br.R);
        X = ps.branch(:,C.br.X);
        y_threshold = 1./(opt.sim.zone1_distance.*abs(R+1j*X));
        relay(:,C.re.threshold) = y_threshold;
        relay(:,C.re.tripped)   = 0;
        relay(:,C.re.id)        = 1:m;
    otherwise
        % make relays of all types
        relay_temp = get_relays(ps,'temperature',opt);
        relay_oc = get_relays(ps,'oc',opt);
        relay_uvls = get_relays(ps,'uvls',opt);
        relay_ufls = get_relays(ps,'ufls',opt);
        relay_dist = get_relays(ps,'distance',opt);
        relay = [relay_temp;
                 relay_oc;
                 relay_uvls;
                 relay_ufls;
                 relay_dist;];
        relay(:,C.re.id)        = 1:(3*m+2*n_shunt);
end

