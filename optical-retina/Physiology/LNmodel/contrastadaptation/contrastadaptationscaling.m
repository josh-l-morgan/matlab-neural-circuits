% contrastadaptationscaling
% this script fit a polynomial to the input-output relation at low contrast
% and than rescales the x-axis of this polynomial to best fit the
% input-output relation measured under high contrast conditions

%% FITTING A CUBIC POLYNOMIAL TO LOW CONTRAST INPUT-OUTPUT
[p,S] = polyfit(lowContrastInput,lowContrastOutput,3);
Y = polyval(p,lowContrastInput);
hold on
plot(lowContrastInput, lowContrastOutput, 'o');
plot(lowContrastInput, Y);

%% RESCALING THE X-AXIS TO FIT HIGH CONTRAST INPUT-OUTPUT

beta = nlinfit(X,y,@xscaling,beta0)
            
        %subfunction here that takes the pvector from polyfit to write a
        %polynomial equation beta should be only an xaxis scaling factor.
