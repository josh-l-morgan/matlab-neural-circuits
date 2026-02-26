function [vec]=crop_vector(vec, element)
%size(vec)

if(length(vec)==1)
	vec=[];
else
	if(element==1)
		vec=vec(2:end);
	elseif(element==length(vec))
		vec=vec(1:end-1);
	else
		vec=[vec(1:element-1) vec(element+1:end)];
	end
end
	