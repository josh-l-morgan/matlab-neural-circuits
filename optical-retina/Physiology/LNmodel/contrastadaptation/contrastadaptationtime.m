xTime = linspace(-500,0,50);
x = linspace(-500,0,500);
scalingOptions = statset ('MaxIter',10000, 'TolFun', 1e-12);
% type = 'ONBSFAST';  %#ok<NASGU>
% type = 'ONBTSLOW';  %#ok<NASGU>
type = 'OFFBSBIPHASIC'; %#ok<NASGU>
% type = 'OFFBTMONOPHASIC';
[m nExperiments] = size(OFFBSBIPHASIC);
%% HIGH CONTRAST
for k=1:nExperiments
    nCells = length(OFFBSBIPHASIC(k).spikes.channel);
    for i=1:nCells
        switch type

            case 'ONBSFAST'
                y = ONBSFAST(k).spikes.highContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [4000, tOne, 10, 2000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                ONBSFAST(k).spikes.highContrastTimeFit(i,:) =...
                    timeFit; %#ok<AGROW>
                trough = x(yhat==min(yhat));
                peak = x(yhat==max(yhat));
                ONBSFAST(k).spikes.highContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>

            case 'ONBTSLOW'
                y = ONBTSLOW(k).spikes.highContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [4000, tOne, 10, 2000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                ONBTSLOW(k).spikes.highContrastTimeFit(i,:) =...
                    timeFit; %#ok<AGROW>
                trough = x(yhat==min(yhat));
                peak = x(yhat==max(yhat));
                ONBTSLOW(k).spikes.highContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>

            case 'OFFBSBIPHASIC'
                y = OFFBSBIPHASIC(k).spikes.highContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [2000, tOne, 10, 4000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                OFFBSBIPHASIC(k).spikes.highContrastTimeFit(i,:) ...
                    = timeFit; %#ok<AGROW>
                trough = x(yhat==max(yhat));
                peak = x(yhat==min(yhat));
                OFFBSBIPHASIC(k).spikes.highContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>

            otherwise
                y = OFFBTMONOPHASIC(k).spikes.highContrastSTAnormalized(i,:);
                tOne = xTime(y==min(y));
                timeFit0 = [-40000, tOne, 15];
                %timeFit = [pOne, tOne, n];
                [timeFit] = nlinfit(xTime, y, @timefieldfit2, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1)));
                OFFBTMONOPHASIC(k).spikes.highContrastTimeFit(i,:) ...
                    = timeFit; %#ok<AGROW>
                deleteFrom = 500 + ceil(x(yhat==min(yhat))) + 1;
                yhat(deleteFrom:end)=[];
                troughLike = x(find(yhat>=0.1*min(yhat),1,'last'));
                peak = x(yhat==min(yhat));
                OFFBTMONOPHASIC(k).spikes.highContrastTimePoint(i,:) = ...
                    [troughLike, peak, peak - troughLike]; %#ok<AGROW>
        end
        clear timeFit* y* peak trough*
    end
end
%% LOW CONTRAST
for k=1:nExperiments
    nCells =length(OFFBSBIPHASIC(k).spikes.channel);
    for i=1:nCells
        switch type

            case 'ONBSFAST'
                y = ONBSFAST(k).spikes.lowContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [4000, tOne, 10, 2000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                ONBSFAST(k).spikes.lowContrastTimeFit(i,:) =...
                    timeFit; %#ok<AGROW>
                trough = x(yhat==min(yhat));
                peak = x(yhat==max(yhat));
                ONBSFAST(k).spikes.lowContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>

            case 'ONBTSLOW'
                y = ONBTSLOW(k).spikes.lowContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [4000, tOne, 10, 2000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                ONBTSLOW(k).spikes.lowContrastTimeFit(i,:) =...
                    timeFit; %#ok<AGROW>
                trough = x(yhat==min(yhat));
                peak = x(yhat==max(yhat));
                ONBTSLOW(k).spikes.lowContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>
                
            case 'OFFBSBIPHASIC'
                y = OFFBSBIPHASIC(k).spikes.lowContrastSTAnormalized(i,:);
                tOne = xTime(y==max(y));
                tTwo = xTime(y==min(y));
                timeFit0 = [2000, tOne, 10, 4000, tTwo];
                %timeFit = [pOne, tOne, n, pTwo, tTwo];
                [timeFit] = nlinfit(xTime, y, @timefieldfit, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1))) -...
                    timeFit(4) * (x./timeFit(5)).^timeFit(3) .*...
                    exp(-timeFit(3)*(x./(timeFit(5)-1)));
                OFFBSBIPHASIC(k).spikes.lowContrastTimeFit(i,:) ...
                    = timeFit; %#ok<AGROW>
                trough = x(yhat==max(yhat));
                peak = x(yhat==min(yhat));
                OFFBSBIPHASIC(k).spikes.lowContrastTimePoint(i,:) = ...
                    [trough, peak, peak - trough]; %#ok<AGROW>

            otherwise
                y = OFFBTMONOPHASIC(k).spikes.lowContrastSTAnormalized(i,:);
                tOne = xTime(y==min(y));
                timeFit0 = [-40000, tOne, 15];
                %timeFit = [pOne, tOne, n];
                [timeFit] = nlinfit(xTime, y, @timefieldfit2, timeFit0,...
                    scalingOptions);
                yhat = timeFit(1) * (x./timeFit(2)).^timeFit(3) .* ...
                    exp(-timeFit(3)*(x./(timeFit(2)-1)));
                OFFBTMONOPHASIC(k).spikes.lowContrastTimeFit(i,:) ...
                    = timeFit; %#ok<AGROW>
                deleteFrom = 500 + ceil(x(yhat==min(yhat))) + 1;
                yhat(deleteFrom:end)=[];
                troughLike = x(find(yhat>=0.1*min(yhat),1,'last'));
                peak = x(yhat==min(yhat));
                OFFBTMONOPHASIC(k).spikes.lowContrastTimePoint(i,:) = ...
                    [troughLike, peak, peak - troughLike]; %#ok<AGROW>
        end
        clear timeFit* y*
    end
end