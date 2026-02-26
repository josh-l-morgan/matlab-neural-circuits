function [filename]=get_temp_filename(prefix, suffix)

counter=1;
filename=[prefix, num2str(counter), suffix];
fstruct=dir(filename);

while(length(fstruct)>0)
	counter=counter+1;
	filename=[prefix, num2str(counter), suffix];
	fstruct=dir(filename);
end	