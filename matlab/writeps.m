function writeps( ps, filename )
% usage: writeps( ps, filename )
% write the power system data file "ps" to "filename"
%  The output is not that pretty, but it is flexible

C  = psconstants;
ps = updateps(ps);

% strip the extension
name = strtok(filename,'.');

% open the file
fid = fopen([filename '.m'],'wt');

%% write the header
fprintf(fid,'function ps = %s\n',name);
fprintf(fid,'\n');

%% write the baseMVA
fprintf(fid,'ps.baseMVA = %f;\n',ps.baseMVA);
fprintf(fid,'ps.casename = %s;\n',ps.casename);
fprintf(fid,'\n');

%% write bus
fprintf(fid,'ps.bus = [...\n');
% column names
fprintf(fid,'%%');
for i=1:length(C.bu.col_names)
    fprintf(fid,'%s ',C.bu.col_names{i});
end
fprintf(fid,'\n');
% data
for i=1:size(ps.bus,1)
	fprintf( fid,' %g',ps.bus(i,:) );
	fprintf( fid,';\n' );
end
fprintf(fid,'];\n\n');

%% write the bus index
%[rows,cols,values] = find(ps.bus_i);
%fprintf(fid,

%% write branch data
fprintf(fid,'ps.branch = [...\n');
% column names
fprintf(fid,'%%');
for i=1:length(C.br.col_names)
    fprintf(fid,'%s ',C.br.col_names{i});
end
fprintf(fid,'\n');
% data
for i = 1:size(ps.branch,1)
	fprintf( fid,' %g',ps.branch(i,:) );
	fprintf( fid,';\n' );
end
fprintf(fid,'];\n\n');

%% write gen
fprintf(fid,'ps.gen = [...\n');
% column names
fprintf(fid,'%%');
for i=1:length(C.ge.col_names)
    fprintf(fid,'%s ',C.ge.col_names{i});
end
fprintf(fid,'\n');
% data
for i = 1:size(ps.gen,1)
	fprintf( fid,' %g',ps.gen(i,:) );
	fprintf( fid,';\n' );
end
fprintf(fid,'];\n\n');

%% write mac
fprintf(fid,'ps.mac = [...\n');
% column names
fprintf(fid,'%%');
for i=1:length(C.ma.col_names)
    fprintf(fid,'%s ',C.ma.col_names{i});
end
fprintf(fid,'\n');
% data
for i = 1:size(ps.mac,1)
	fprintf( fid,' %g',ps.mac(i,:) );
	fprintf( fid,';\n' );
end
fprintf(fid,'];\n\n');

%% write shunt
if isfield(ps,'shunt')
    fprintf(fid,'ps.shunt = [...\n');
    % column names
    fprintf(fid,'%%');
    for i=1:length(C.sh.col_names)
        fprintf(fid,'%s ',C.sh.col_names{i});
    end
    fprintf(fid,'\n');
    % data
    for i = 1:size(ps.shunt,1)
        fprintf( fid,' %g',ps.shunt(i,:) );
        fprintf( fid,';\n' );
    end
    fprintf(fid,'];\n\n');
end

%% areas
if isfield(ps,'areas')
    fprintf(fid,'ps.areas = [...\n');
    % data
    for i = 1:size(ps.areas,1)
        fprintf( fid,' %g',ps.areas(i,:) );
        fprintf( fid,';\n' );
    end
    fprintf(fid,'];\n\n');
end

%% gencost
if isfield(ps,'gencost')
    fprintf(fid,'ps.gencost = [...\n');
    % data
    for i = 1:size(ps.gencost,1)
        fprintf( fid,' %g',ps.gencost(i,:) );
        fprintf( fid,';\n' );
    end
    fprintf(fid,'];\n\n');
end


