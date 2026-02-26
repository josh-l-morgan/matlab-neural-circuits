function spikes = contrastadaptationparameters(spikes, i)

%highContrastLateParameters = [sigma, mu, maximumFiring]
scalingOptions = statset ('MaxIter', 10000);
muSeed = find(spikes(i).highContrastLateOutput <= ...
    (max(spikes(i).highContrastLateOutput) - min(spikes(i).highContrastLateOutput))/2, 1, 'last');
highContrastLateParameters0 = [0.1, spikes(i).highContrastLateInput(muSeed), max(spikes(i).highContrastLateOutput)];
[highContrastLateParameters] = nlinfit(spikes(i).highContrastLateInput,...
    spikes(i).highContrastLateOutput, @contrastadaptationfitone, highContrastLateParameters0, scalingOptions);
spikes(i).highContrastLateParameters = highContrastLateParameters;

%highContrastEarlyParameters = [sensitivity maintainedDrive]
passInput.sigma = highContrastLateParameters(1);
passInput.mu = highContrastLateParameters(2);
passInput.maximumFiring = highContrastLateParameters(3);
passInput.data = spikes(i).highContrastEarlyInput;

highContrastEarlyParameters0 = [1,0];
[highContrastEarlyParameters] = nlinfit(passInput,...
    spikes(i).highContrastEarlyOutput, @contrastadaptationfittwo, highContrastEarlyParameters0, scalingOptions);
spikes(i).highContrastEarlyParameters = highContrastEarlyParameters;


%lowContrastEarlyParameters = [sensitivity maintainedDrive]
passInput.data = spikes(i).lowContrastEarlyInput;

lowContrastEarlyParameters0 = [1,0];
[lowContrastEarlyParameters] = nlinfit(passInput,...
    spikes(i).lowContrastEarlyOutput, @contrastadaptationfittwo, lowContrastEarlyParameters0, scalingOptions);
spikes(i).lowContrastEarlyParameters = lowContrastEarlyParameters;

%lowContrastLateParameters = [sensitivity maintainedDrive]
passInput.data = spikes(i).lowContrastLateInput;

lowContrastLateParameters0 = [1,0];
[lowContrastLateParameters] = nlinfit(passInput,...
    spikes(i).lowContrastLateOutput, @contrastadaptationfittwo, lowContrastLateParameters0, scalingOptions);
spikes(i).lowContrastLateParameters = lowContrastLateParameters;








