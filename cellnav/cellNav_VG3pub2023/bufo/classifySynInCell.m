function[synIds] = classifySynInCell(syn,s,cid);


%%Collect synapses
for i = 1:length(s)
    hit = [1:length(syn.pre)]';

    if ~isempty(s(i).synType)
        hit = intersect(hit,find(syn.synType == s(i).synType));
    else
        hit = intersect(hit,find(syn.synType<inf)) ;
    end

    if isempty(s(i).input)
        part = syn.pre(hit);

    elseif s(i).input
        hit = intersect(hit, find(syn.post == cid));
        part = syn.pre(hit);

    else
        hit = intersect(hit, find(syn.pre == cid));
        part = syn.post(hit);
    end


    if ~s(i).checkCids
        synIds{i} = hit;
    else
        isPart = part * 0;
        for h = 1:length(hit)
            if sum(s(i).cids == part(h))
                isPart(h) = 1;
            end
        end
        synIds{i} = hit(isPart>0);
    end
    
end




