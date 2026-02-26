function[sProp] = maxPredProp(pathL,prop,L)

%% spread properties along pred looking at L in two directions

if ~exist('L','var');
    L = 5;
end

sProp = prop;
for t = 1:length(pathL.tips)
    
    node = pathL.tips(t);
    clear predProp predNode meanProp
    predProp = prop(node);
    predNode = node;
    for p = 1:length(pathL.pred)
        if node==pathL.bases(t)
            break
        end
        predProp(p) = prop(node);
        predNode(p) = node;
        
        node = pathL.pred(node);
        if node==0
            break
        end
    end
    
    meanProp = zeros(length(predProp),1);
    for p = 1 : length(predProp)
        grab = [p-L:p+L];
        grab = grab((grab>=1) & (grab<=length(predProp)));
%         if length(grab)==0
%             return
%         end
        meanProp(p) =  max(predProp(grab));
        sProp(predNode(p)) = meanProp(p);
        grabNum(predNode(p)) = length(grab);
    end
    %sProp(predNode) = t;
end