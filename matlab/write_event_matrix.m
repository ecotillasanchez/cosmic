function write_event_matrix(event,filename)
% writes out an event matrix for later usage

% strip the extension
name = strtok(filename,'.');

% open the file
fid = fopen([name '.m'],'wt');

% get constants
C  = psconstants;

%% write the header
fprintf(fid,'function event = %s\n',name);
fprintf(fid,'\n');

%% write the matrix
fprintf(fid,'event = [...\n');
% column names
fprintf(fid,'%%');
for i=1:length(C.ev.col_names)
    fprintf(fid,'%s ',C.ev.col_names{i});
end
fprintf(fid,'\n');
% data
for i=1:size(event,1)
	fprintf( fid,' %g',event(i,:) );
	fprintf( fid,';\n' );
end
fprintf(fid,'];\n\n');
