function[links] = breakLoop(links,uPairs,jointVol)

% allNodes = unique(links);
% nonPairs = setdiff(uPairs,allNodes);
% uPairs = cat(1,uPairs,nonPairs);
[uPairs ix] = sort(uPairs);
jointVol = jointVol(ix);

[paths, showP] = getPaths(links);

%% find loop
startOver = 1;
while startOver
    startOver = 0;
    for p = 2:length(paths)  %run all paths
        for l = 1:p-1   %compare to previous paths
            if startOver == 0
                shared = intersect(paths{l},paths{p});
                shared = shared(shared>0);
                if length(shared)>1
                    %get paths in loop
                    path1 = paths{p};
                    pRange = find((path1==shared(1)) | (path1 == shared(2)));
                    path1 = path1(pRange(1):pRange(end));
                    path2 = paths{l};
                    pRange = find((path2==shared(1)) | (path2 == shared(2)));
                    path2 = path2(pRange(1):pRange(end));

                    %Get weakes links
                    segs1  =  [path1(1:end-1) path1(2:end)];
                    segs2 = [path2(1:end-1) path2(2:end)];
                    segs = cat(1,segs1,segs2);
                    nodes = sort(unique(segs));
                    strength = jointVol(ismember(uPairs, nodes));
                    weakest = find(strength == min(strength),1);
                    rm = nodes(weakest); %node to remove
                    [n s] = find(segs == rm);
                    spouse = segs(sub2ind(size(segs),n,~(s-1)+1));
                    weakLinks = [ones(length(spouse),1)*rm  spouse];

                    for r = 1:size(weakLinks,1) %remove weakest
                        links = links(~((links(:,1) == rm) & (links(:,2) == spouse(r))),:);
                        links = links(~((links(:,2) == rm) & (links(:,1) == spouse(r))),:);
                    end

                    %start over
                    [paths, showP] = getPaths(links);
                    startOver = 1;
                    break
                end
            end
        end
    end

end



