function spikes = contrastadaptationparameters(spikes, i)

%highContrastParameters = [sigma, mu, maximumFiring]
scalingOptions = statset ('MaxIter', 10000);
muSeed = find(spikes(i).highContrastOutput <= ...
    (max(spikes(i).highContrastOutput) - min(spikes(i).highContrastOutput))/2, 1, 'last');
highContrastParameters0 = [0.1, spikes(i).highContrastInput(muSeed), max(spikes(i).highContrastOutput)];
[highContrastParameters] = nlinfit(spikes(i).highContrastInput,...
    spikes(i).highContrastOutput, @contrastadaptationfitone, highContrastParameters0, scalingOptions);
spikes(i).highContrastParameters = highContrastParameters;

%highContrastEarlyParameters = [sensitivity maintainedDrive]
passInput.sigma = highContrastParameters(1);
passInput.mu = highContrastParameters(2);
passInput.maximumFiring = highContrastParameters(3);
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








