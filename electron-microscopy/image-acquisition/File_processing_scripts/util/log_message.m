% EMNET: Parallel implementation of convolutional networks
% Developed and maintained by Viren Jain <viren@mit.edu>
% Do not distribute without permission.

function [fid]=log_message(n, message)
if(isempty(n))
	n.params.save_string='/home/viren/emnet';
end

message=['[',datestr(now,0),'] ',message];
fprintf(1, [message,'\n']);
fid=fopen([n.params.save_string, 'log'], 'a+');
for i=1:size(message,1)
	fprintf(fid, [message(i,:), '\n']);
end
fclose(fid);

