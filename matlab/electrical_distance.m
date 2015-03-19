function E = electrical_distance(Ybus,method,use_cluster,use_pf_jacobian,V)
% calculate a power flow electrical distance matrix
% usage: E = electrical_distance(Ybus,method,use_cluster,use_pf_jacobian,V)
% 
% methods available:
%  1. Jacobian: | d(Ti-Tj)/d(Pi-Pj) |, where Ti is the bus angle in radians
%  2. Ybus:     | inv(Ybus) |
%  3. 'Efficient' uses the correct electrical distance method, based on dP_dTheta
% the 'use_cluster' flagg will write a set of n scripts, in which the columns
% of E1 are calculated (not E = s

n = size(Ybus,1);
% default inputs:
if nargin<2, method='Efficient'; end
if nargin<3, use_cluster=false; end
if nargin<4, use_pf_jacobian=false; end
if nargin<5, V = ones(n,1); end

% constants:
j = 1i;

% prep work
if use_cluster
    work_dir = pwd;
    out_dir = [work_dir sprintf('/E_matrix_%d',n)];
    mkdir(out_dir);
    % copy the m file to out_dir
    system(sprintf('cp diag_of_inverse.m %s',out_dir));
end

% initialize the output
E = zeros(n,n);

% pre-calculate the Jacobian
if ~use_pf_jacobian
    Ibus = Ybus * V;
    diagV     = spdiags(V, 0, n, n);
    diagIbus  = spdiags(Ibus, 0, n, n);
    dSbus_dVa = j * diagV * conj(diagIbus - Ybus * diagV);
    dP_dTheta = real(dSbus_dVa);
else
    [~,~,dP_dTheta,~] = fnPowerFlowJac(Ybus,real(V),imag(V));
end

% do the remaining calculations
switch method
    case 'FullJac'
        
    case 'Efficient' % an efficient version of the Jacobian method
        G = dP_dTheta;
        [mn,r] = min(abs(angle(V)));
        if mn>1e-12
            error('Reference bus should be zero angle');
        end
        not_r = true(n,1);
        not_r(r) = false;
        Ginv_sub = inv(G(not_r,not_r));
        Ginv_diag = diag(Ginv_sub);
        E(r,not_r) = Ginv_diag';
        E(not_r,r) = Ginv_diag;
        M = ones(n-1,1)*Ginv_diag';
        E(not_r,not_r) = M + M' - Ginv_sub - Ginv_sub';
        
    case 'Jacobian' % Jacobian method
        disp('Calculating electrical distance...');
        disp('This is a stupid way to do this...');
        
        % calculate the distances
        
        if ~use_cluster % simple version of the code:
            E1 = zeros(n,n);
            for b = 1:n
                tic;
                % calculate the buses in question
                bus_set = setdiff(1:n,b);
                % re-calculate the jacobian
                J = dP_dTheta(bus_set,bus_set);
                % Find the diagonal of the inverse of J
                E1(bus_set,b) = diag_of_inverse(J);
                if n>1000
                    fprintf('Row %d of %d took %5.2f sec\n',b,n,toc);
                end
            end
            % Calculate the distance matrix
            E = sqrt(E1.*E1');
        else % linux cluster version of the code:
            % split the job into n_workers sub-jobs
            n_workers = 45;
            cols_per_worker = ceil(n/n_workers);
            wd = pwd;
            cd(out_dir);
            % calculate the matrix for each row/column:
            disp('Writing jacobians to files');
            nJ = 0;
            for b=1:n
                J_fname = sprintf('J_%05d.mat',b);
                if ~exist(J_fname,'file')
                    % calculate the buses in question
		    J = dP_dTheta;
		    J(:,b) = [];
		    J(b,:) = [];
                    % save it to a file
                    save(J_fname,'J');
                    nJ = nJ+1;
                end
            end
            fprintf('Wrote %d jacobians\n',nJ);
            last_col = 0;
            for worker = 1:n_workers
                fprintf('Preparing job for worker #%d\n',worker);
                first_col = last_col + 1;
                last_col  = min(first_col + cols_per_worker,n);
                if last_col<first_col
                    break;
                end
                % produce a worker job to calculate the columns
                qsub_diag_of_inverse(first_col:last_col);
                % print something
                fprintf('Submitted job for rows %d through %d\n',first_col,last_col);
            end
            cd(wd);
            E = out_dir;
        end
    case 'invYbus' % inv(Ybus) method
        if use_cluster, error('use_cluster doesn''t work with Ybus method'); end
        E = abs(inv(Ybus));
    case 'Impedance'
        if use_cluster, error('use_cluster doesn''t work with impedance method'); end
        
        % do the LU decomposition on Ybus
        disp('Calculating Zbus...');
        if n<1000
            Zbus = full(inv(Ybus));
        else
            [L,U,p,q,R] = lu(Ybus,'vector');
            Zbus = zeros(n,n);
            E = zeros(n,n);
            tic;
            for i=1:n
                % if it has been more than 10 seconds, tell the user what is going on.
                if toc>10
                    fprintf(' column %5d of %5d\n',i,n);
                    tic;
                end
                % calculate one column of the Zbus matrix
                x = sparse(i,1,1,n,1);
                Zbus(:,i) = lusolve(x,L,U,p,q,R);
            end
        end
        disp('Calculating distance...');
        tic;
        for a=1:n
            % if it has been more than 10 seconds, tell the user what is going on.
            if toc>10 
                fprintf(' row %5d of %5d\n',a,n);
                tic;
            end
            all_but_a = setdiff((1:n),a);
            for b=all_but_a
                Zab = Zbus([a b],[a b]);
                det_Zab = det(Zab);
                if abs(det_Zab)<1e-12
                    E(a,b) = abs(imag(Zab(1,2)));
                    keyboard
                else
                    y_ab = -Zab(1,2)/det_Zab;
                    %E(a,b) = abs(1/y_ab);
                    E(a,b) = abs(imag(1/y_ab));
                    E(a,b) = abs(1/imag(y_ab));
                end
            end
        end
        E = sqrt(E.*E');
    otherwise
        error('Unknown method');
end

if use_cluster
end

%% functions...

function qsub_diag_of_inverse(cols)
% create and execute a script for submission via qsub

%% prep work:
jobname = sprintf('run_%05d_%05d',cols(1),cols(end));

%% build the m-file for this job
disp('Writing the m files');
m_fname = sprintf('%s.m',jobname);
mf = fopen(m_fname,'w');
if isempty(mf)
    error('Could not create m file');
end
n_Ei = 0;
for i=cols
    Ei_filename = sprintf('Ei_%05d.mat',i);
    if ~exist(Ei_filename,'file')
        fprintf(mf,sprintf('if ~exist(''%s'',''file'') \n',Ei_filename));
        fprintf(mf,sprintf('  load J_%05d; \n',i));
        fprintf(mf,sprintf('  Ei = diag_of_inverse(J); \n'));
        fprintf(mf,sprintf('  save %s Ei; \n',Ei_filename));
        fprintf(mf,sprintf('  disp(''Finished calculating column %d'');\n',i));
        fprintf(mf,sprintf('end\n'));
        n_Ei = n_Ei+1;
    end
end
fclose(mf);

if n_Ei==0
    fprintf('All of columns %d to %d have been calculated\n',cols(1),cols(end));
    return
end

%% launch the job
runtime_minutes = 5*length(cols);
runtime_hours = runtime_minutes/60;
memory_gb = 4;
vmemory_gb = 5;
run_on_vacc(m_fname,true,runtime_hours,memory_gb,vmemory_gb);


