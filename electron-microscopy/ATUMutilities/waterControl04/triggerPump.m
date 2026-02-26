function[] = triggerPump

global watProps


%% Make sound
if watProps.soundOn
    y = sin([1:10000*watProps.pumpDuration]/10);
    sound(y,10000)
end


%% open serial
%seriallist

try
    if strcmp(watProps.ps.Status,'closed')
        watProps.ps.Terminator='CR';
        fopen(watProps.ps);
    elseif strcmp(watProps.ps.Status,'open')
        watProps.ps.Terminator='CR';
    end
catch err
    watProps.ps = serial('COM3');
    watProps.ps.Terminator='CR';
    fopen(watProps.ps);
end

pumpcmd='RUN';

%% dummy
'pumping it'
%disp(sprintf('pump duration %0.02f ',watProps.pumpDuration))

for rep = 1:watProps.pumpDuration
    fprintf(watProps.ps,pumpcmd);
    if watProps.pumpDuration>1
        pause(3);
    end
end

%fclose(watProps.ps)


