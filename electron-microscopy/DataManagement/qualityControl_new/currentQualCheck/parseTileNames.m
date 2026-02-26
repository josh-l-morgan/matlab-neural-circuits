function[mos] = parseTileNames(picNams)

%% Find digits and prefixes
for i = 1:length(picNams)
    nam = picNams{i};
    %% Generic parser
    dig = regexp(nam,'\d');
    start = 1; c = 0;
    lastDig = -100;
    vals = []; pre = {};
    for d = 1:length(dig)
        if dig(d)~= (lastDig + 1)
            c = c + 1;
            val = nam(dig(d));
            pre{c} = nam(start:dig(d)-1);
        else
            val = [val nam(dig(d))];
        end
        vals(c) = str2num(val);
        start = dig(d) + 1;
        lastDig = dig(d);
    end
    pnam(i).pre = pre;
    pnam(i).vals = vals;
end

%% Match Patterns
pat = 1:length(pnam);
for i = 1:length(pnam)
    if pat(i) ==i;
        pre = pnam(i).pre;
        for s = i+1:length(pnam)
            toMatch = pnam(s).pre;
            same = 0;
            if length(toMatch) == length(pre)
                for p = 1:length(pre)
                    if strcmp(pre{p},toMatch{p}) 
                        same = same + 1;
                    end %if same pre
                end %check all pres
            end %if same length
            if same == length(pre)
                pat(s) = pat(i);
            end
        end %check subsequent names
    end %if pat not already assigned
end%check all names

%% get dominant pattern
upats = unique(pat);
for p = 1:length(upats)
    sumPats(p) = sum(pat == upats(p));
end
usePats = find(pat==upats(find(sumPats == max(sumPats),1)));

%%  assign possitions

template = pnam(usePats(1)).pre;
for t = 1:length(template)
    for u = 1:length(usePats)
        tVals(u,t) = pnam(usePats(u)).vals(t);
    end
        numVals(t) = length(unique(tVals(:,t)));  
end

valList = sort(numVals,'descend');
useVals = find(numVals>= valList(2));
useVals = useVals(1:2);

if regexp(template{useVals(1)},'[cx]','ignorecase')
    useVals = fliplr(useVals);
elseif  regexp(template{useVals(2)},'[ry]','ignorecase')
    useVals = fliplr(useVals);
end

mos.Pos = [tVals(:,useVals(1))  tVals(:,useVals(2))];
mos.Dims = max(mos.Pos,[],1);
mos.Files = usePats;
end





