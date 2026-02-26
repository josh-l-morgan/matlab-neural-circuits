function[I] = sumSub(sub,dim,fsize);

if ~exist('fsize','var')
    fsize = max(sub,[],1);
    fmin = min(sub,[],1);
    fsize = fsize-fmin+1;
else
    fmin = [0 0 0];
end

if ~exist('dim','var')
    dim = 1;
end

if dim == 1
    dims = [3 2];
elseif dim == 2
    dims = [3 1];
elseif dim == 3
    dims = [1 2];
end

            I = zeros(fsize(dims));
                inds = sub2ind(fsize(dims),sub(:,dims(1))-fmin(dims(1))+1,sub(:,dims(2))-fmin(dims(2))+1);
                uinds = unique(inds);
                
                if length(uinds)>1
                    hinds = hist(inds,uinds);
                else
                    hinds = length(inds);
                end
                
                I(uinds) = I(uinds) + hinds';
                % newHeight(uinds) = sub(:,dim);
                
                
                