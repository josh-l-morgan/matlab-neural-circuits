function[synV] = getSynapsePropertyForType(sms)

%%Generate number of rgcsubtypes or bipolar cell subtypes innervated by 
%%VG3s in a way that reduces the influence of local clustering

global tis

clear synV

if ~exist('minDist','var')
  minDist = 0;
end  
if ~exist('checkTypeID','var')
    checkTypeID = 1;
end

%count synapses between pairs of VG3s and RGCs with minDist exclusion length
synStack = [];
for v = 1:length(sms)
    sm = sms(v).sm;
    syn = sm.syn;
    synStack = cat(1,synStack,[syn.pre syn.post zeros(length(syn.pre),1) + sm.cid sm.syn2Skel.closest ]);
end

%% Get types

isID = (synStack(:,1)>0) & (synStack(:,2)>0);
synStack = synStack(isID,:);

clear typeStack subTypeStack
for i = 1:size(synStack,1)
    cid2 = synStack(i,2);
    targ = find(tis.cids==cid2);
    if isempty(targ)
        typeStack(i) = 0;
        subTypeStack(i) =    0;
    else
        typeStack(i) = tis.cells.type.typeID(targ);
        subTypeStack(i) =   tis.cells.type.subTypeID(targ);
    end
end


synV.pre = synStack(:,1);
synV.post = synStack(:,2);
synV.vCid = synStack(:,3);
synV.synID = synStack(:,4);
synV.type = typeStack(:);
synV.subTypeStack = subTypeStack(:);
for i = 1:length(typeStack)
    if subTypeStack(i);
        synV.subTypeNames{i} = tis.cells.type.subTypeNames{typeStack(i)}{subTypeStack(i)};
    else
        synV.subTypeNames{i} = 'none';
    end
end





