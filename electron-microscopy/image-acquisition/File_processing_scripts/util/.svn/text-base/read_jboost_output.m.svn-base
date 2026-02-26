function [set_info]=read_jboost_output(fname)

fid=fopen(fname, 'r');
text=textscan(fid, '%n%n%n%n%n%n', 'delimiter', ':','Headerlines',1,'EndOfLine','\n');
set_info=zeros([size(text{1},1) 2]);

for i=1:size(set_info,1)
	set_info(i,1)=text{2}(i);
	set_info(i,2)=text{3}(i);
	set_info(i,3)=text{4}(i);
end