

wif = GetMyWafer;
focSec = checkFocus(wif);
%focSec = manFocus(wif,focSec);
retake = shapeOutput(wif,focSec);
retakeControl(wif,retake);