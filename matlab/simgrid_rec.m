function [outputs,ps] = simgrid_rec(ps,event,outfilename,opt)
% usage: [outputs,ps] = simgrid_rec(ps,event,outfilename,opt)
% perform a cascading failure simulation (in rectangular formulation)
% inputs:
%  ps           - powersystem structure. see psconstants.
%  event        - a matrix defining the exogenous events to simulate. see psconstants.
%  outfilename  - the name for the output file
%  opt          - an option structure. see psoptions.
% outputs:
%  outputs.         - a structure of output data, including the following:
%    .success       - a bool flag indicating whether we were able to simulate to the end.
%    .t_simulated   - vector of time periods in the final results.
%    .outfilename   - the name of the file where the results are stored.
%    .demand_lost   - the amount of demand that was shed during the simulation.
%    .computer_time - the amount of wall-clock time required to perform this simulation
%    .endo_events   - matrix that logs endogenous events during the simulation.
%  ps               - the powersystem structure at the end of the simulation.
 

%% process the inputs
if nargin<2 || isempty(event)
    if ~isempty(ps.event) 
        event = ps.event; 
    else
        error('simgrid:err','undefined event matrix.');
    end
    if ~isfield('branch_i')
        error('simgrid:err','need ps.branch_i to run simgrid');
    end
end
if nargin<3 || isempty(outfilename)
    outfilename = 'outfile.csv';
end
if nargin<4 || isempty(opt)
    opt = psoptions;
end

%% clean up the ps structure
ps = updateps(ps);

%% edit the output file name
start_time = datestr(now,30);
outfilename = sprintf('%s_%s.csv',strtok(outfilename,'.'),start_time);
tracefilename = sprintf('trace_%s_%s.mat',strtok(outfilename,'.'),start_time);

%% prepare the outputs
outputs.success         = false;
outputs.t_simulated     = [];
outputs.outfilename     = outfilename;
outputs.demand_lost     = 0;
outputs.computer_time   = [];   ct = clock;
outputs.start_time      = start_time;     
ps.event_record         = [];
event_record            = [];

%% other prep work
C           = psconstants;
verbose     = opt.verbose;
discrete    = false;

%% check the event matrix
event = sortrows(event);
if event(1,C.ev.type)~=C.ev.start
    error('simgrid:err','first event should be a start event.');
end
t_0   = event(1,1);
event = event(2:end,:); % remove the start event
if event(end,C.ev.type)~=C.ev.finish
    error('simgrid:err','last event should be a finish event.');
end
t_end = event(end,1);
t = t_0;

%% print something
out = fopen(outfilename,'w');
if isempty(out)
    error('simgrid:err','could not open outfile: %s',outfilename);
end
fprintf(out,'starting simulation at t = %g\n',t);
fclose(out);

if verbose
    fprintf('starting simulation at t = %g\n',t);
    fprintf('writing results to %s\n',outfilename);
end

%% get Ybus, x, y
% build the Ybus
if ~isfield(ps,'Ybus') || ~isfield(ps,'Yf') || ~isfield(ps,'Yt')
    [ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
end
% build x and y
[x,y] = get_xy_rec(ps,opt);

%% step through the simulation
event_no = 1;
while t < t_end
    % process exogenous discrete events at this time period
    while event(event_no,1) <= t
        [ps,new_discrete] = process_event(ps,event(event_no,:),opt);
        event_no = event_no+1;
        discrete = discrete | new_discrete;     % for simultaneous events
    end
    if discrete
        t = t + opt.sim.t_eps;
    end
    % find the next time interval
    t_next = event(event_no,1);
    % try to simulate between t and t_next
    if opt.verbose
        fprintf('\n Simulating from t=%g s to t=%g s\n',t,t_next);
    end
    [ps,t_out,X,Y] = simgrid_interval_rec(ps,t,t_next,x,y,opt);
    if opt.verbose
        fprintf(' Completed simulation up until t=%g\n',t_out(end));
    end
    
    % log to files
    if opt.sim.writelog
        if opt.verbose
            fprintf(' Writing simulation results to %s\n',outfilename);
        end
        write_state(outfilename,t_out,X,Y);
    end;
    
    % update time and solutions for next interval
    t = t_next;
    x = X(:,end);
    y = Y(:,end);
    % record events
    event_record            = [event_record; ps.event_record]; %#ok<AGROW>
    ps.event_record         = [];
end

%% update outputs and clean up
if t_out(end) < (t_end - 1)
    outputs.success         = false;
    outputs.demand_lost     = sum(ps.shunt(:,C.sh.P));
else
    outputs.success         = true;
    outputs.demand_lost     = sum(ps.shunt(:,C.sh.P)) - ps.shunt(:,C.sh.P)'*ps.shunt(:,C.sh.factor);
end
outputs.t_simulated     = [t_0 t_out(end)];
if opt.sim.writelog
    outputs.outfilename     = outfilename;
else 
    outputs.outfilename     = [];
end
outputs.event_record    = event_record;
outputs.computer_time   = etime(clock,ct);
save(tracefilename,'x','y');
if opt.verbose
    fprintf('Completed simulation from %d sec. to %d sec. \n',t_0,t_next);
end
