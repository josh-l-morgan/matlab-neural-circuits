
function[difNum reportString] = checkNewNetwork(synMat,makeSyn,touchMat);


checkTouch = (makeSyn>0) - (touchMat>0);
badSyn = sum(checkTouch(:)>0);

mistakeNum = sum(sum(synMat,2) ~= sum(makeSyn,2));
badSpine = sum(sum(makeSyn,1) ~= sum(synMat,1));

synDif = sum(synMat(:)) - sum(makeSyn(:));

reportString = sprintf('%d non touch synapses,%d difference in synapse number, %d axon differences, %d spine differences', badSyn,synDif,mistakeNum,badSpine);

difNum = [badSyn synDif mistakeNum badSpine];