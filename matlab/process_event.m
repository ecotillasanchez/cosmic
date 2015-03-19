function [ps,discrete] = process_event(ps,event,opt)
% usage: [ps,discrete] = process_event(ps,event,opt)
% process the event specified in the event vector (one row of an event matrix)


if size(event,1)>1
    error('process_event:err','can only process one event at a time.');
end

verbose = opt.verbose;
discrete = false;   % is this a discrete event?

% extract some data
C = psconstants;
t = event(C.ev.time);

% record the event
ps.event_record = [ps.event_record;event];

% do the processing
switch event(C.ev.type)
    case C.ev.start
        if verbose, fprintf(' simulation start event.'); end
        
    case C.ev.fault     % three phase to ground fault by default
        error('process_event:err','faults not supported yet.');
        
    case C.ev.trip_branch
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 0
            ps.branch(branch_ix,C.br.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d tripped...\n',t,branch_id); end
        end
        
    case C.ev.close_branch
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 1
            ps.branch(branch_ix,C.br.status) = 1;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d closed...\n',t,branch_id); end
        end
        
    case C.ev.trip_bus
        bus_no = event(1,C.ev.bus_loc);
        br_set = any( ps.branch(:,1:2)==bus_no, 2 );
        ps.branch(br_set,C.br.status) = 0;
        % trip gens and shunts at this bus
        ps.gen(ps.gen(:,1)==bus_no,C.gen.status) = 0;
        ps.shunt(ps.shunt(:,1)==bus_no,C.shunt.status) = 0;
        discrete = true;
        if verbose, fprintf('  t = %.4f: Bus %d tripped...\n',t,bus_no); end
        
    case C.ev.trip_gen
        gen_id = event(C.ev.gen_loc);
        gen_ix = ps.gen_i(gen_id);
        if ps.gen(gen_ix,C.gen.status)~=0
            ps.gen(gen_ix,C.gen.status) = 0;
            ps.gen(gen_ix,C.Pg) = 0;
            ps.gen(gen_ix,C.Qg) = 0;
            if verbose, fprintf('  t = %.4f: Gen %d tripped...\n',t,gen_id); end
            discrete = true;
        end
        
    case C.ev.shed_load
        shunt_id = event(C.ev.shunt_loc);
        change_by = event(C.ev.change_by);
        prev_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
        if change_by
            ps.shunt(shunt_id,C.sh.factor) = max(ps.shunt(shunt_id,C.sh.factor) - event(C.ev.quantity), 0);
            curr_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        else
            curr_load = max(prev_load - event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.P) = max(ps.shunt(shunt_id,C.sh.P)-event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.Q) = max(ps.shunt(shunt_id,C.sh.Q)-event(C.ev.quantity),0);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        end
        
    case C.ev.trip_shunt
        shunt_id = event(C.ev.shunt_loc);
        bus_no = ps.shunt(shunt_id,1);
        if ps.shunt(shunt_id,C.sh.status) ~= 0
            ps.shunt(shunt_id,C.sh.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: shunt %d at bus %d tripped...\n',t,shunt_id,bus_no); end
        end
        
    case C.ev.close_shunt
        shunt_id = event(C.ev.shunt_loc);
        shunt_id = ps.shunt_i(shunt_id);
        bus_no = ps.shunt(shunt_id,1);
        if ps.shunt(shunt_id,C.sh.status) ~= 1
            ps.shunt(shunt_id,C.sh.status) = 1;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d at bus %d closed...\n',t,shunt_id,bus_no); end
        end
    case C.ev.oc_relay
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 0
            ps.branch(branch_ix,C.br.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d tripped...\n',t,branch_id); end
        end
    case C.ev.temp_relay
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 0
            ps.branch(branch_ix,C.br.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d tripped...\n',t,branch_id); end
        end
    case C.ev.dist_relay
        branch_id = event(C.ev.branch_loc);
        branch_ix = ps.branch_i(branch_id);
        if ps.branch(branch_ix,C.br.status) ~= 0
            ps.branch(branch_ix,C.br.status) = 0;
            discrete = true;
            % output the branch id
            if verbose, fprintf('  t = %.4f: Branch %d tripped...\n',t,branch_id); end
        end
    case C.ev.uvls_relay
        shunt_id = event(C.ev.shunt_loc);
        change_by = event(C.ev.change_by);
        prev_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
        if change_by
            ps.shunt(shunt_id,C.sh.factor) = max(ps.shunt(shunt_id,C.sh.factor) - event(C.ev.quantity), 0);
            curr_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        else
            curr_load = max(prev_load - event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.P) = max(ps.shunt(shunt_id,C.sh.P)-event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.Q) = max(ps.shunt(shunt_id,C.sh.Q)-event(C.ev.quantity),0);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        end
    case C.ev.ufls_relay
        shunt_id = event(C.ev.shunt_loc);
        change_by = event(C.ev.change_by);
        prev_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
        if change_by
            ps.shunt(shunt_id,C.sh.factor) = max(ps.shunt(shunt_id,C.sh.factor) - event(C.ev.quantity), 0);
            curr_load = ps.shunt(shunt_id,C.sh.P).*ps.shunt(shunt_id,C.sh.factor);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        else
            curr_load = max(prev_load - event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.P) = max(ps.shunt(shunt_id,C.sh.P)-event(C.ev.quantity),0);
            ps.shunt(shunt_id,C.sh.Q) = max(ps.shunt(shunt_id,C.sh.Q)-event(C.ev.quantity),0);
            discrete = true;
            bus_no = ps.shunt(shunt_id,1);
            if opt.verbose, fprintf('  t = %.4f: %.2f MW of load shedding at bus %d...\n',t,prev_load-curr_load,bus_no); end
        end
        
        
    otherwise
        error('process_event:err','Unknown event type');
end
