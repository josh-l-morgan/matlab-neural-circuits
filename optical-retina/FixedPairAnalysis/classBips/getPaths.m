function[paths,showP] = getPaths(links);


%% specifically check for loops

%%Check everything
check = ones(size(links,1),1);
paths = {};
p = 1;
while sum(check)  %check all possible links
    targ = find(check,1);
    target = links(targ,1);  %node to check
    paths{p} = target;
    %%Fill path one
    while p <= length(paths);  %while there are more paths to explore
        while 1; % if still finding spouses
            path = paths{p};
            target = path(end);
            [seg side]  = find(links == target);
            ok = check(seg)>0;
            oside = ~(side-1)+1;
            spouse = links(sub2ind(size(links),seg(ok),oside(ok)));
            check(seg) = 0;
            
            if sum(ok) == 0, break,end  %finished path if no one new
            if sum(ok) == 1,
             paths{p} = cat(1,paths{p},spouse(1)); %add to current path
            else
                for np =1:length(spouse)  %seed new paths
                    paths{length(paths)+1} = [target; spouse(np)];
                end  %end add new spouses
            end
        end % if still finding spouses
        p = p + 1;
    end %while there are more paths

end   % check all possible nodes


%% show paths
if nargout>1
clear showP
for i = 1: length(paths)
    path = paths{i};
    showP(1:length(path),i) = path;
end
end