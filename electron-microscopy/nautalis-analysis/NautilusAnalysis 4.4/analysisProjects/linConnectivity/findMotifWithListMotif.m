
if 0
    clear all
    load('MPN.mat');
    load([MPN 'obI.mat']);
    targCell = 125;
    maxDist = 5;
    iptsetpref('ImshowBorder','tight');
    
    load('D:\LGNs1\Analysis\sm.mat')
    
    
end

%% analyze Motifs


% save(['D:\LGNs1\Analysis\motifImages\listMotifWorkSpaceR.mat'])
% load(['D:\LGNs1\Analysis\motifImages\listMotifWorkSpaceR.mat'])

isTarg =(sm.pre == targCell) | (sm.post == targCell);
isRGC = sm.preClass == 1;
isTC = sm.postClass == 2;
isPreLIN = sm.preClass == 3;
isPostLIN = sm.postClass == 3;
isPreUNK = sm.preClass == 4;


%search all RGC input, check for non RGC TC
%search all LIN input, check for isolated LIN output
%search for all UNK input


%% find checkIn relative motifs



[checkRgc] = checkMotifAroundSyns(sm,find(isRGC & isTarg),maxDist);
[checkPreLin] = checkMotifAroundSyns(sm,find(isPreLIN & (sm.post == targCell)),maxDist);
[checkUnk] = checkMotifAroundSyns(sm,find(isPreUNK & isTarg),maxDist);
[checkPostLin] = checkMotifAroundSyns(sm,find(isPostLIN & (sm.pre == targCell)),maxDist);
[checkTc] = checkMotifAroundSyns(sm,find(isTC & isTarg),maxDist);

sum([checkRgc.m.syn2skelDi]>0)
sum([checkRgc.m.diadTCs]>0)
sum([checkRgc.m.preTcDi]>0)



motifCodes = [...
    '1 = RGC only, 2 = TC only, 3 = LINin only, 4 = LINout only, 5 = UNKin only, ' ...
    '6 = recip, 7 = autapse,'...
    '8 = RGC-LIN-TC diad, 9 = RGC-LIN-TC triad, 10 = RGC-LIN-TC mult triad, ' ...
    '11 = RGC-LIN-LIN diad, 12 = RGC-LIN-LIN IN triad,  13 = RGC-LIN-LIN out triad, 14 = RGC-LIN-LIN recip triad, ' ...
    '15 = LIN1-LIN-TC diad (top), 16 = LIN-LIN1-TC diad (mid), '...
    '17 = LIN1-LIN-TC triad(top), 18 = LIN-LIN1-TC triad(mid), 19 = LIN-LIN-TC recip triad, '...
    '20 = LLL diad, 21 = LLL tri top, 22 = LLL tri mid, 23 = LLL tri bot, 24 = LLL recip triad, '...
    '25 = UNK-LIN-TC(diad), 26 = UNK-LIN-TC triad, 27 = LIN-TC without RGC, 28 = LIN-TC with near RGC, diad'];



pos = [];
for i = 1:length(checkPostLin.m)
    m = checkPostLin.m(i)
    if m.goodPos
        if m.reciprocal
            pos = [pos;checkPostLin.pos(i,:)];
        end
    end
end

anchors = pos;
dSamp =  (obI.em.res .* [4 4 1])./1000;
%dSamp =  (obI.em.res .* [4 4 1])./1000./obI.em.dsRes;
%dSamp = dSamp ./ [4 4 2];



anchors(:,1) = anchors(:,1)/dSamp(1);
anchors(:,2) = anchors(:,2)/dSamp(2);
anchors(:,3) = anchors(:,3)/dSamp(3);
%anchors = round(anchors);
anchors(anchors<1) = 1;

%%
for i = 1:size(anchors,1)
    i
    sprintf('%.0f %.0f %.0f',anchors(i,2),anchors(i,1),anchors(i,3))
    pause 
    
end


