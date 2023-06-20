function [t_steps,X,Y,Z] = solve_dae(f,g,h,aux,x0,y0,t_span,opt)
% usage: [t_steps,X,Y,Z] = solve_dae(f,g,h,aux,x0,y0,t_span,opt)
% this function uses the first order trapezoidal rule to solve
% the dae defined by the following:
%
%  f - the differential variables
%  g - the algebraic variables
%  h - the integer variables
%  aux - the glocal id for relays or ix depends on the inputs
%  x0 - the initial set of differential variables
%  y0 - the initial set of algebraic variables
%  t_span - time span
%  opt - options

% process default inputs and outputs
if nargin<8, opt = numerics_options; end;
if isempty(h), check_discrete=false; else check_discrete=true; end
Z = [];

global t_delay t_prev_check dist2threshold 

% constants, sizes and defaults
max_newton_its  = opt.sim.max_iters; % maximum number of newton iterations
abs_tol         = opt.sim.tolerance; % newton convergence tolerance
dif_tol_min     = 0.01;              % increase step size
dif_tol_max     = 0.05;              % decrease step size
dt_min          = 0.005;             % default min step size
dt_max          = 1;                 % default max step size
alpha_0         = 1;
min_alpha       = 2^-10;
nx              = length(x0);
X               = x0;
Y               = y0;
xy0 = [x0;y0];

unloop_delta    = 0;                 % 1 = unloop the delta so that it is within [-2*pi,2*pi]
opt.sim.uvls_tdelay_ini = 0.5;       % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.5;       % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;       % 1 sec delay for dist relay.
opt.sim.emergency_control = 0;
% find starting point
t0 = t_span(1);

% default time step
t_final = t_span(end);
dt0     = dt_max*0.10;
t_steps = t0;
% get the initial state of the dae system
f0 = f(t0,x0,y0);
if check_discrete
    z0 = h(t0,xy0,[]);
    z0_prev = z0; % initialize z0_prev as z0
end

% do the LU on the jacobian to get the symbolic work done
% [L,U,p,q,R] = lu(jacobian);

% set initial step size
dt = dt0;

while t0<t_final
    % choose a guess for f,g at next point
    t1 = t0 + dt;
    dt_next = [];
    if t1 > t_final
        t1 = t_final;
        dt = t1-t0;
    end  % integrate up to t_final
    x1 = x0 + f0*dt;    % could use second order estimate here
    y1 = y0;            % could solve for something more intelligent
    % print something
%     fprintf(' t = %f sec. and dt = %f sec...\n',t0,dt);

    % iteratively update x1 and y1 until it converges
    newton_it = 0;
    while newton_it < max_newton_its+1
        % get the next residuals, jacobian blocks
        [f1,df_dx1,df_dy1] = f(t1,x1,y1);
        [g1,dg_dx1,dg_dy1] = g(t1,x1,y1);
        % check the mismatch and quit
        trapz_mismatch  = [x0 - x1 + (dt/2).*f0 + (dt/2).*f1;   % f portion
            g1];                                                % g portion
        max_trap_mis = max(abs(trapz_mismatch));
%         fprintf('   Mismatch = %g\n',max_trap_mis);
        if max_trap_mis < abs_tol
            %fprintf('.converged in %d iterations with alpha=%g ...\n',newton_it,alpha);
            % if it converged, first check for any events at t0<th<t1
            % [outputs] = get_event_thresholds(inputs)
            break;
        end
        % bail out if we are at max number of iterations
        if newton_it == max_newton_its
%             break
            error(' time step did not converge...');
        end
        
        % build the jacobian that will be used in Newton's method
        trapz_df_dx = -speye(nx) + (dt/2).*df_dx1;
        trapz_df_dy = (dt/2).*df_dy1;
        jacobian    = [trapz_df_dx trapz_df_dy;
            dg_dx1 dg_dy1];
        
        % find the search direction
        p = -(jacobian\trapz_mismatch);
        % implement the backstepping to get the Newton step size
        alpha = alpha_0;
        while 1
            x1 = x1 + alpha.*p(1:nx);
            y1 = y1 + alpha.*p((nx+1):end);
            f1 = f(t1,x1,y1);
            f_mis = x0 - x1 + (dt/2).*f0 + (dt/2).*f1;
            g1 = g(t1,x1,y1);
            old_trap_mis = max_trap_mis;
            max_trap_mis = max(max(abs(f_mis)),max(abs(g1)));
            
            if max_trap_mis < old_trap_mis
%                 fprintf('Newton step completed.\n');
                break
            else
                fprintf('Reducing Newton step size.\n');
                alpha = alpha/2;
            end
            if alpha<min_alpha
                fprintf('Algorithm failure in the newton step.\n');
                return
            end
        end
        newton_it = newton_it + 1;
    end
    
    % at this point assume that the solution is a good one.
    solution_good = true;

    dt_adjust_relay = false;
    % check if a threshold was crossed
    if check_discrete
        xy1 = [x1;y1];
        z1 = h(t1,xy1,[]);
        % check to see if we hit a threshold
        ix              = aux([],1);
        hit_thresh = abs(z1)<opt.sim.eps_thresh;
        down_crossed = z1 <= -opt.sim.eps_thresh;
        up_crossed = z1(1:ix.re.nrelay) >= opt.sim.eps_thresh; % exclude emgency control events
        % figure out when the threshold was crossed
        % adjust dt for down_crossed
        if ~any(hit_thresh) && any(down_crossed)
            if any(sign(z1(down_crossed))~=sign(z0_prev(down_crossed)))
                crossed = find(down_crossed);
                crossing = crossed(sign(z1(down_crossed))~=sign(z0_prev(down_crossed)));
                dz = z1(crossing) - z0_prev(crossing);
                t_thresh = -(dt./dz).*z0_prev(crossing) + t0;
                [~,which_thresh] = min(t_thresh);
                dt = max(t_thresh(which_thresh) - t0,dt_min);
                solution_good = false;
                dt_adjust_relay = true;
            end
        end
        % adjust dt for up_crossed
        if ~any(hit_thresh) && any(up_crossed)
            if any(sign(z1(up_crossed))~=sign(z0_prev(up_crossed)))
                crossed = find(up_crossed);
                crossing = crossed(sign(z1(up_crossed))~=sign(z0_prev(up_crossed)));
                dz = z1(crossing) - z0_prev(crossing);
                t_thresh = -(dt./dz).*z0_prev(crossing) + t0;
                [~,which_thresh] = min(t_thresh);
                dt = max(t_thresh(which_thresh) - t0,dt_min);
                solution_good = false;
                dt_adjust_relay = true;
            end
        end
        % adjust dt in case dt is greater than the remaining t_delay
        if any(down_crossed(1:ix.re.nrelay))
            crossed = down_crossed(1:ix.re.nrelay);
            t_remain = t_delay(aux(crossed,0));
            if dt > max(min(t_remain),dt_min/5)
                dt = max(min(t_remain),dt_min/5);
                solution_good = false;
                dt_adjust_relay = true;
            end
        end
        z0_prev = z0;
        z0 = z1;
    end
        % adjust step size if needed
    if opt.sim.var_step
        if (max(abs(f1-f0)) >= dif_tol_max) && (dt > dt_min)
            % reduce step size and go back (ie trash x1, f1, etc)
            dt = max( dt/2, dt_min);
            solution_good = false;
            dt_next = min(2*dt, dt_max);
        else
            if (max(abs(f1-f0)) < dif_tol_min) && ~dt_adjust_relay
                % increase next step size
                dt_next = min(2*dt, dt_max);
            end
        end
    end
    % record the solution for this step
    if solution_good
        % commit the results to memory
        % if it converged and no events, then save t1, x1 and y1; advance
        t_prev      = t0;
        t0          = t1;           % commit time advance
        [z1,~,~,Temperature] = h(t1,xy1,dt);          % compute dist2threshold and update state_a
        
        if ~isempty(dt_next)
            dt          = dt_next;
        end
        % unloop delta signals within [-2pi,2pi]
        if unloop_delta
            ix              = aux([],1);  %#ok<*UNRCH>
            x1(ix.x.delta)  = rem(x1(ix.x.delta),2*pi);
        else
            ix              = aux([],1);
        end
        
        % record the temperature signal in X variable
        x1(ix.x.temp)  = Temperature; % in X, temperature is referenced to ambient temperature 20 degree C
        
        x0          = x1;
        y0          = y1;
        f0          = f1;
        % save the variables
        X           = [X x1];       %#ok<*AGROW>
        Y           = [Y y1];
        t_steps     = [t_steps;t1];

        % when X and Y are committed, update the time delay for each relay
            down_crossed_re = down_crossed;
            hit_thresh_re   = hit_thresh;
            
        if any(hit_thresh_re) || any(down_crossed_re)
            % take down the time when first time hit or crossed
            relay_event = or(hit_thresh_re,down_crossed_re);
            [rows,~] = find(relay_event);
            n_relay = size(rows,1);
            for i = 1:n_relay
                % if the trace has never acrossed its threshold or its t_prev_check has
                % been restored, t_prev_check is NaN
                if isnan(t_prev_check(aux(rows(i),0)))
                    t_prev_check(aux(rows(i),0)) = t_prev;
                end
            end
            % when it is below threshold, reduce the time delay until it hits 0
            if any(sign(z1(down_crossed_re)) == sign(z0_prev(down_crossed_re)))
                crossed = find(down_crossed_re);
                stay_crossed = crossed(sign(z1(down_crossed_re))==sign(z0_prev(down_crossed_re))); % find the local id for the 'down_crossed' signals
                n_stay_crossed = size(stay_crossed,1);
                for i = 1:n_stay_crossed
                    if ismember(stay_crossed(i),[ix.re.oc])
                        t_delay(aux(stay_crossed(i),0)) = dist2threshold(aux(stay_crossed(i),0))/(-z1(stay_crossed(i)));
                    else
                        int = t1 - t_prev_check(aux(stay_crossed(i),0));
                        t_delay(aux(stay_crossed(i),0)) = t_delay(aux(stay_crossed(i),0)) - int; % reduce its time delay according to its global id
                        t_prev_check(aux(stay_crossed(i),0)) = t1;
                    end
                end
                % when t_delay = 0, break and process this event
                if any(t_delay(aux(stay_crossed,0))<=0)
                    tripped = stay_crossed(t_delay(aux(stay_crossed,0))<=0);
                    Z = false(size(relay_event,1),1);
                    Z(tripped) = true;
                    t_delay(t_delay<0)=0;
                    break
                end
            end
        end
        % when it is above threshold, increase the time delay until it hits
        % full range of its setting.
        if any(up_crossed)
            if any(sign(z1(up_crossed)) == sign(z0_prev(up_crossed)))
                crossed = find(up_crossed);
                stay_crossed = crossed(sign(z1(up_crossed))==sign(z0_prev(up_crossed)));    % find the local id for the 'up_crossed' signals
                n_stay_crossed = size(stay_crossed,1);
                for j = 1:n_stay_crossed
                    if ~ismember(stay_crossed(j),[ix.re.oc])
                        if ~isnan(t_prev_check(aux(stay_crossed(j),0)))
                            int = t1 - t_prev_check(aux(stay_crossed(j),0));
                            t_delay(aux(stay_crossed(j),0)) = t_delay(aux(stay_crossed(j),0)) + int; % incease its time delay according to its global id
                            t_prev_check(aux(stay_crossed(j),0)) = t1;
                            % when t_delay restores to initial value, set the
                            % previous check time back to NaN
                            if ismember(stay_crossed(j),[ix.re.uvls])
                                if t_delay(aux(stay_crossed(j),0)) >= opt.sim.uvls_tdelay_ini
                                    t_delay(aux(stay_crossed(j),0)) = opt.sim.uvls_tdelay_ini;
                                    t_prev_check(aux(stay_crossed(j),0)) = nan;
                                end
                            elseif ismember(stay_crossed(j),[ix.re.ufls])
                                if t_delay(aux(stay_crossed(j),0)) >= opt.sim.ufls_tdelay_ini
                                    t_delay(aux(stay_crossed(j),0)) = opt.sim.ufls_tdelay_ini;
                                    t_prev_check(aux(stay_crossed(j),0)) = nan;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
