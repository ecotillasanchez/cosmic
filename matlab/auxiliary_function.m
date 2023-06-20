function [ output ] = auxiliary_function(local_id,ix,ps,opt)
% usage: [ global_relay_id ] = auxiliary_funciton(local_id,ix,ps,opt)
% output is the global relay id if the inputs are (local_id, 0)
% output is the ix if the inputs are ([], 1)

C = psconstants;

if ~ix && ~isempty(local_id)
    output = ps.relay(local_id,C.re.id);
    
elseif ix && isempty(local_id)
    n           = size(ps.bus,1);
    ng          = size(ps.mac,1);
    m           = size(ps.branch,1);
    n_sh        = size(ps.shunt,1);
    output      = get_indices(n,ng,m,n_sh,opt);
end
end