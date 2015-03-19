function ps_large = replicate_case(ps_small,S)
% usage: ps_large = replicate_case(ps_small,S)
% replicates a ps case within a matrix in order to obtain a larger ps case
C           = psconstants;
ps_small    = updateps(ps_small);

if nargin<2
    S = 10; % square matrix size (e.g. S=10 would convert case39 into case3900)
end

M = ones(S,S);  
alpha = 100;    % renumbering scale (it has to be larger than the number of buses in ps_small)

% dummy branch to interconnect matrix elements (FIXME, hard-coded for a generic high voltage branch)
dummy.br = [0 0 0.0005 0.0272 0 9900 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0];

% interconnecting buses (FIXME, hard-coded for case39)
north_west = 2;     north_east = 29;
south_west = 9;     south_east = 23;

% iterate through the matrix and compose ps_large, always connect left/up
for col = 1:S
    for row = 1:S
        if col == 1 && row == 1
            % do not connect anything, initialize ps_large
            ps      = ps_small;
            r_this  = alpha;
            % renumber stuff
            ps.bus(:,1)         = r_this+ps.bus(:,1);
            ps.branch(:,1:2)    = r_this+ps.branch(:,1:2);
            ps.gen(:,1)         = r_this+ps.gen(:,1);
            ps.shunt(:,1)       = r_this+ps.shunt(:,1);
            ps.mac(:,1)         = r_this+ps.mac(:,1);
            ps.exc(:,1)         = r_this+ps.exc(:,1);
            ps.gov(:,1)         = r_this+ps.gov(:,1);
            ps_large = ps;
        else
            ps                        = ps_small;
            ps.bus(31,C.bus.type)     = C.PV;           % unique slack
            ps.gen(2,C.gen.type)      = C.PV;
            ps.gen(2,C.gen.part_fact) = 0;
            r_this  = alpha*sub2ind(size(M),row,col);   % bus renum idx
            % renumber stuff 
            ps.bus(:,1)         = r_this+ps.bus(:,1);
            ps.branch(:,1:2)    = r_this+ps.branch(:,1:2);
            ps.gen(:,1)         = r_this+ps.gen(:,1);
            ps.shunt(:,1)       = r_this+ps.shunt(:,1);
            ps.mac(:,1)         = r_this+ps.mac(:,1);
            ps.exc(:,1)         = r_this+ps.exc(:,1);
            ps.gov(:,1)         = r_this+ps.gov(:,1);
            % append to ps_large
            ps_large.bus    = [ps_large.bus     ;ps.bus];
            ps_large.branch = [ps_large.branch  ;ps.branch];
            ps_large.gen    = [ps_large.gen     ;ps.gen];
            ps_large.shunt  = [ps_large.shunt   ;ps.shunt];
            ps_large.mac    = [ps_large.mac     ;ps.mac];
            ps_large.exc    = [ps_large.exc     ;ps.exc];
            ps_large.gov    = [ps_large.gov     ;ps.gov];
            % now interconnect in the appropiate directions
            if col == 1 && row ~= 1
                % connect just up
                r_up    = alpha*sub2ind(size(M),row-1,col);  % up bus renum idx
                % connect dummy branches
                % north_west
                tie_this = ps_large.bus(:,1) == r_this+north_west;
                tie_up   = ps_large.bus(:,1) == r_up+south_west;
                dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_up,1)];
                ps_large.branch = [ps_large.branch  ;dummy.br];
                % north_east
                tie_this = ps_large.bus(:,1) == r_this+north_east;
                tie_up   = ps_large.bus(:,1) == r_up+south_east;
                dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_up,1)];
                ps_large.branch = [ps_large.branch  ;dummy.br];   
            else if col ~= 1 && row == 1
                    % connect just left
                    r_left  = alpha*sub2ind(size(M),row,col-1);  % up bus renum idx
                    % connect dummy branches
                    % north_west
                    tie_this = ps_large.bus(:,1) == r_this+north_west;
                    tie_left = ps_large.bus(:,1) == r_left+north_east;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_left,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                    % south_west
                    tie_this = ps_large.bus(:,1) == r_this+south_west;
                    tie_left = ps_large.bus(:,1) == r_left+south_east;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_left,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                else
                    % connect up and left
                    r_up    = alpha*sub2ind(size(M),row-1,col);  % up bus renum idx
                    r_left  = alpha*sub2ind(size(M),row,col-1);  % up bus renum idx
                    % connect dummy branches
                    % north_west
                    tie_this = ps_large.bus(:,1) == r_this+north_west;
                    tie_up   = ps_large.bus(:,1) == r_up+south_west;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_up,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                    % north_east
                    tie_this = ps_large.bus(:,1) == r_this+north_east;
                    tie_up   = ps_large.bus(:,1) == r_up+south_east;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_up,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                    % north_west
                    tie_this = ps_large.bus(:,1) == r_this+north_west;
                    tie_left = ps_large.bus(:,1) == r_left+north_east;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_left,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                    % south_west
                    tie_this = ps_large.bus(:,1) == r_this+south_west;
                    tie_left = ps_large.bus(:,1) == r_left+south_east;
                    dummy.br(1,1:2) = [ps_large.bus(tie_this,1) ps_large.bus(tie_left,1)];
                    ps_large.branch = [ps_large.branch  ;dummy.br];
                end
            end
        end
    end
end
ps_large.areas   = [];
ps_large.gencost = [];

% reset device id's
ps_large.branch(:,C.br.id)  = 0;
ps_large.gen(:,C.ge.id)     = 0;
ps_large.shunt(:,C.sh.id)   = 0;
ps_large = updateps(ps_large);
