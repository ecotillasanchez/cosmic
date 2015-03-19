function [relay_type,relay_location,relay_index] = get_relay_event_info(relay_event,ps)
% usage: [relay_type,relay_location] = get_relay_event_info(relay_event,ps)
% gets the relay event information (type and location)
C = psconstants;

relay_type = ps.relay(relay_event,C.re.type);
relay_id = find(relay_event);
num_relay = size(relay_id,1);
relay_index = zeros(num_relay,1);
relay_location = zeros(num_relay,1);

for i = 1:num_relay
    if relay_type(i) == C.relay.temp
        relay_ind = relay_id(i);
        location = ps.relay(relay_ind,C.relay.branch_loc);
    elseif relay_type(i) == C.relay.oc
        relay_ind = relay_id(i);
        location = ps.relay(relay_ind,C.relay.branch_loc);
    elseif relay_type(i) == C.relay.dist
        relay_ind = relay_id(i);
        location = ps.relay(relay_ind,C.relay.branch_loc);
    elseif relay_type(i) == C.relay.uvls
        relay_ind = relay_id(i);
        location = ps.relay(relay_ind,C.relay.shunt_loc);
    elseif relay_type(i) == C.relay.ufls
        relay_ind = relay_id(i);
        location = ps.relay(relay_ind,C.relay.shunt_loc);
    end
    relay_location(i) = location;
    relay_index(i) = relay_ind;
end


