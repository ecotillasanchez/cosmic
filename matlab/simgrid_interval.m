function [ps,t_out,X,Y] = simgrid_interval(ps,t,t_next,x0,y0,opt)
% usage: [ps,t_out,X,Y] = simgrid_interval(ps,t,t_next,x0,y0,opt)
% integrate PS DAEs from t to t_next, recursively if endogenous events
%
% inputs:
%  ps - power system structure, see psconstants.
%  t  - initial simulation time.
%  t_next - end simulation time assuming no relay events.
%  x0 and y0 - current state of the system
%  opt - options inherited from simgrid
%
% outputs:
%  ps - the power systems structure at the end of the simulation.
%   .endo_events   - matrix that logs endogenous events during the simulation.
%  t_out - vector with times where the system was evaluated at
%  X and Y - state of the system during the integration period

% constants and data
C           = psconstants;
n           = size(ps.bus,1);
n_macs      = size(ps.mac,1);
m           = size(ps.branch,1);
n_shunts    = size(ps.shunt,1);
j = 1i;

angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

if COI_weight
    weight     = ps.mac(:,C.ma.M);
else
    weight     = ps.gen(:,C.ge.mBase);
end

% check whether the last relay event splitted the system
br_status           = ps.branch(:,C.br.status) == C.CLOSED;
is_first_subgraph   = findSubGraphs(ps.bus(:,C.bu.id), ps.branch(br_status,C.br.f:C.br.t),1);

if ~all(is_first_subgraph)
    % the network partitioned
    if opt.verbose
        fprintf('  t = %.4f: The network partitioned into two islands...\n',t);
    end
    % get ps structures and DAE variables for the two subnets
    net1 = logical(is_first_subgraph); net2 = ~net1;
    ps1  = subsetps(ps,net1);   [x01,y01] = get_xy(ps1,opt);
    ps2  = subsetps(ps,net2);   [x02,y02] = get_xy(ps2,opt);

    % step down a recursion level and solve for the two subnets
    [ps1,t_out1,X1,Y1] = simgrid_interval(ps1,t,t_next,x01,y01,opt);
    [ps2,t_out2,X2,Y2] = simgrid_interval(ps2,t,t_next,x02,y02,opt);

    % aggregate the outputs from the lower recursion levels
    [ps,t_out,X,Y] = superset_odeout(ps,ps1,ps2,t_out1,t_out2,X1,Y1,X2,Y2,opt);

elseif (n_macs == 0 || n_shunts == 0 || (ps.shunt(:,C.sh.P)'*ps.shunt(:,C.sh.factor) == 0) || size(ps.bus(:,1),1) == 1)
    % there is no generation or load in this network, or it is just one disconnected bus
    t_out           = t;
    X               = nan(size(x0));
    if isempty(X)                           % disconnected bus without mac
        X = [X; nan(1,size(X,2))];
    end
    Y               = nan(size(y0));
    if ~isempty(ps.shunt), ps.shunt(:,C.sh.factor) = 0; end
else
    % nothing special (for now), try to integrate the network DAEs from t to t_next
    ix          = get_indices(n,n_macs,m,n_shunts,opt);
    mac_bus_i   = ps.bus_i(ps.mac(:,1));
    temp_ref    = 0;

    % create a temp reference bus if the subgrid doesn't have one
    if ~any(ps.bus(:,C.bu.type) == 3)
        [~,ref] = max(ps.gen(:,C.ge.Pmax));
        temp_ref = ps.bus(ismember(ps.bus(:,1),ps.gen(ref,1)),C.bu.type);
        ps.bus(ismember(ps.bus(:,1),ps.gen(ref,1)),C.bu.type) = C.REF;
    end
    
    % updata the frequency signal sources for every load bus
    ps = update_load_freq_source(ps);
    
    % recalculate algebraic variables according to the updated Ybus
    [ps.Ybus,ps.Yf,ps.Yt,ps.Yft,ps.sh_ft] = getYbus(ps,false);    
    y_new = solve_algebraic(t,x0,y0,ps,opt);
    
    if isempty(y_new)
        % could not solve for the algebraic variables, shut down subgrid
        t_out           = t;
        X               = nan(size(x0));
        if isempty(X) % disconnected bus without mac
            X = [X; nan(1,size(X,2))];
        end
        Y               = nan(size(y0));
        if ~isempty(ps.shunt), ps.shunt(:,C.sh.factor) = 0; end
        % remove temp reference bus, if needed
        if temp_ref, ps.bus(ismember(ps.bus(:,1),ps.gen(ref,1)),C.bu.type) = temp_ref; end
        
    else
        % we should be able to start integrating the DAE for this subgrid
        xy0 = [x0;y_new];
        y0 = y_new;
        
        % choose integration scheme
        switch opt.sim.integration_scheme
            case 1
                % trapezoidal rule
                clear fn_f; clear fn_g; clear fn_h;
                fn_f = @(t,x,y) differential_eqs(t,x,y,ps,opt);
                fn_g = @(t,x,y) algebraic_eqs(t,x,y,ps,opt);
                fn_h = @(t,xy,dt) endo_event(t,xy,ix,ps,dt,opt);
                fn_aux= @(local_id,ix) auxiliary_function(local_id,ix,ps,opt);
                [t_ode,X,Y,Z] = solve_dae(fn_f,fn_g,fn_h,fn_aux,x0,y0,t:opt.sim.dt_default:t_next,opt);
                XY_ode = [X;Y]';
                
            case 2
                % implicit, ode15i
                clear fn_fg;
                xyp0 = zeros(size(xy0));
                event_handle = @(t,xy,xyp) endo_event_i(t,xy,xyp,ix,ps);
                fn_fg = @(t,xy,xyp) differential_algebraic_i(t,xy,xyp,ps,opt);
                options = odeset(  'Jacobian', @(t,xy,xyp) get_jacobian_i(t,xy,xyp,ps,opt), ...
                    'Events', event_handle, ...
                    'Stats','off');
                [xy0,xyp0]                      = decic(fn_fg,t,xy0,zeros(size(xy0)),xyp0,zeros(size(xyp0)));
                [t_ode,XY_ode]                  = ode15i(fn_fg,t:opt.sim.dt_default:t_next,xy0,xyp0,options);
                
            case 3
                % explicit, ode15s
                clear fn_fg;
                event_handle = @(t,xy) endo_event(t,xy,ix,ps);
                fn_fg = @(t,xy) differential_algebraic(t,xy,ix.nx,ps,opt);
                mass_matrix = sparse(1:ix.nx,1:ix.nx,1,ix.nx+ix.ny,ix.nx+ix.ny);
                options = odeset(   'Mass',mass_matrix, ...
                    'MassSingular','yes', ...
                    'Jacobian', @(t,xy) get_jacobian(t,xy,ix.nx,ps,opt), ...
                    'MaxOrder', 2, ...
                    'BDF', 'off', ...
                    'Stats','off', ...
                    'RelTol', 1e-3, ...
                    'AbsTol', 1e-6, ...
                    'Events', event_handle, ...
                    'NormControl','off');
                [t_ode,XY_ode]                  = ode15s(fn_fg,t:opt.sim.dt_default:t_next,xy0,options);
        end
        
        % organize data from the DAE output
        t_out       = t_ode';
        x_rows      = (1:ix.nx);
        y_rows      = (ix.nx+1):(ix.nx+ix.ny);
        X           = XY_ode(:,x_rows)';
        Y           = XY_ode(:,y_rows)';
        x_end       = X(:,end);
        y_end       = Y(:,end);
        
        % store DAE data back into the ps structure
        mac_Thetas  = y_end(ix.y.theta(mac_bus_i));
        deltas      = x_end(ix.x.delta);
        if ~angle_ref
            delta_sys   = y_end(ix.y.delta_sys);
            delta_m     = deltas + delta_sys - mac_Thetas;
        else
            delta_coi   = sum(weight.*deltas,1)/sum(weight,1);
            delta_m     = deltas - delta_coi - mac_Thetas;   % Revisit this
        end
        
        ps.mac(:,C.mac.delta_m)     = delta_m;
        ps.mac(:,C.mac.omega)       = x_end(ix.x.omega_pu)*2*pi*ps.frequency;
        ps.mac(:,C.mac.Pm)          = x_end(ix.x.Pm);
        ps.mac(:,C.mac.Eap)         = x_end(ix.x.Eap);
        ps.exc(:,C.ex.E1)       	= x_end(ix.x.E1);
        ps.exc(:,C.ex.Efd)      	= x_end(ix.x.Efd);
        ps.gov(:,C.go.P3)          	= x_end(ix.x.P3);
        ps.bus(:,C.bus.Vmag)        = y_end(ix.y.Vmag);
        ps.bus(:,C.bus.Vang)        = y_end(ix.y.theta)*180/pi;
        ps.relay(ix.re.temp,C.re.state_a) = x_end(ix.x.temp);
        
        % branch results
        Vmag = y_end(ix.y.Vmag);
        Theta = y_end(ix.y.theta);
        V = Vmag.*exp(j.*Theta);
        If = ps.Yf*V; % branch status is accounted for in Yf
        It = ps.Yt*V; % branch status is accounted for in Yt
        F = ps.bus_i(ps.branch(:,1));
        T = ps.bus_i(ps.branch(:,2));
        Sf = V(F) .* conj(If);
        St = V(T) .* conj(It);
        ps.branch(:,C.br.Imag_f) = abs(If);
        ps.branch(:,C.br.Imag_t) = abs(It);
        ps.branch(:,C.br.Pf) = real(Sf) * ps.baseMVA;
        ps.branch(:,C.br.Qf) = imag(Sf) * ps.baseMVA;
        ps.branch(:,C.br.Pt) = real(St) * ps.baseMVA;
        ps.branch(:,C.br.Qt) = imag(St) * ps.baseMVA;        

        % branch power loss
        Yft = - ps.Yft;
        sh_ft = ps.sh_ft;
        Sft = V(F).*conj((V(F)-V(T)).*Yft + V(F).*sh_ft);
        Stf = V(T).*conj((V(T)-V(F)).*Yft + V(T).*sh_ft);
        lineloss = Sft + Stf;
        negvalue = real(lineloss)<0;
        lineloss(negvalue) = lineloss(negvalue) * -1;
        ps.branch(:,C.br.lineloss) = lineloss*ps.baseMVA;
        
        % preparation work for the emergency control
        % calculate resulting ZIPE load after the powerflow
        % get the load bus injections with a ZIPE model
        D               = ps.bus_i(ps.shunt(:,1));
        S_load_base     = (ps.shunt(:,C.sh.P) + j*ps.shunt(:,C.sh.Q)).*ps.shunt(:,C.sh.factor)/ps.baseMVA;
        S_load_P        = S_load_base.*ps.shunt(:,C.sh.frac_S);
        Sd              = sparse(D,3,S_load_P,n,5);
        S_load_Z        = S_load_base.*ps.shunt(:,C.sh.frac_Z);
        Sd              = Sd + sparse(D,1,S_load_Z,n,5);
        S_load_I        = S_load_base.*(1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E)));
        Sd              = Sd + sparse(D,2,S_load_I,n,5);
        S_load_E        = S_load_base.*ps.shunt(:,C.sh.frac_E);
        Sd              = Sd + sparse(D,4,S_load_E,n,5);
        S_load_E_gamma  = ps.shunt(:,C.sh.gamma);
        Sd              = Sd + sparse(D,5,S_load_E_gamma,n,5);
        zipe_cols       = 5;   % assuming it is ZIPE model for now
        if zipe_cols == 1
            S_zipe = Sd;
        elseif zipe_cols == 5
            S_Z = Sd(:,1) .* Vmag.^2;
            S_I = Sd(:,2) .* Vmag;
            S_P = Sd(:,3);
            S_E = Sd(:,4) .* Vmag.^Sd(:,5);
            S_zipe = S_Z + S_I + S_P + S_E;
        else
            error('zipe load model matrix is not the right size');
        end
        
        ps.shunt(:,C.sh.current_P) = real(S_zipe(ps.shunt(:,C.sh.type)==1))* ps.baseMVA ;
        ps.shunt(:,C.sh.current_Q) = imag(S_zipe(ps.shunt(:,C.sh.type)==1))* ps.baseMVA ;
        
        % remove temp reference bus, if needed
        if temp_ref, ps.bus(ismember(ps.bus(:,1),ps.gen(ref,1)),C.bu.type) = temp_ref; end
            
            relay_event = Z;
            % if there was a relay event, apply it.
            if ~isempty(relay_event)
                
                t_event = t_out(end);
                % process reley event
                [ps] = process_relay_event(t_event,relay_event,ps,opt);
                t_down = t_event + opt.sim.t_eps;
                % step down a recursive level to continue solving from the event to t_next
                [ps,t_out_down,X_down,Y_down] = simgrid_interval(ps,t_down,t_next,x_end,y_end,opt);
                
                % merge the outputs from this recursion level and the one immediately below
                t_out       = [t_out t_out_down];
                X           = [X X_down];
                Y           = [Y Y_down];
                
                % reset relay_event as []
                relay_event = []; %#ok<NASGU>
            end
            
    end
end

return
