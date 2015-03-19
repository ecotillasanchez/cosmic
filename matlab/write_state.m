function write_state(outfilename,t_vector,X,Y)
% usage: write_state(outfilename,t_vector,X,Y)
% writes the state variables (t,x,y) to the output file

fid = fopen(outfilename,'a+');
if isempty(fid), error('write_state:err','could not open %s .',outfilename); end

for k = 1:length(t_vector)
    fprintf(fid,'%.8f,',t_vector(k));
    fprintf(fid,'%g,',X(:,k));
    fprintf(fid,'%g,',Y(:,k));
    fprintf(fid,'\n');
end
fclose(fid);
