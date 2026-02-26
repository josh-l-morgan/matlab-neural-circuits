








%{

Core relationships:
length constant is (square  root of (membrane resistance /axial
resistance)), sqrt(rm/ra)

specific membrane resistence (Rm) is the inverse of conductance, Ohms *
cm^2, the inverse of conductance

membrane resistance (rm) = Rm / (2 * pi * r)

axial resistance (ra = Ri / (pi * r^2), Ri is the specific resistance of
cytoplasm in oms * cm. 


Rm = 1/Cm;
ra = Ri / (pi * r^2)
rm = Rm / (2 * pi * r)
lc = sqrt(rm/ra)



NEURON default parameters:
Ri = 100 = axial resistance in Ohm * cm
Cm = 1 = membrane capacitance in uf/cm^2
Gm = 5 * 10^-4 = membrane conductance in S/cm^2

Conductances to test:
5 * 10 ^4, 1 * 10^4, 0.5 * 10^4




An axonal length constant can be
calculated using the mean biophysical parameters (s mem cap = 0.01 pF * um^-2; 
internal g = 1.3 mSuS*um^-1; mem g = 0.83 pS*um^-2), a constant diameter of 1.6 um,
and the simplifying assumption of a semi-infinite cable. The resulting value is 791 um,


%}

%%
r = 0.5;
logGm = [-4 -2.75 -2.5 -2];
Gm = 10.^ logGm;
Cm = 1;
Ri = 100;

rCm = r / 10^4;
Rm = 1./Gm;
ra = Ri / (pi * rCm^2);
rm = Rm ./ (2 * pi * rCm);
lcCm = sqrt(rm./ra);


%lcCm = sqrt((2*rCm*Rm)/(4 * Ri));

lc = lcCm * 10^4
clf
plot(logGm,lc)






%%
Gm = 5*10^4; %test specific conductance for membrane S/cm^2
r = .5; %test radius for neurite in micrometers
rCm = r / 10^4; %convert micrometers to centimeters
Ri = 100; %test axial resistance Ohm * cm


Rm = 1 / Gm; %Calculate specific membrane resistance from specific membrane conductance
rm = Rm / (2 * pi * rCm); % Calculate membrane resistance from specific membrane resitance and radius
ra = Ri / (pi * rCm^2); % Calculate axial resistance from specific axial resistance and radius
lcCM = sqrt(rm/ra); % Calculate length constant (cm) from membrane resistance and axial resistance

lc = lcCM * 10^4 % Convert centimeters to micrometers




Gm = 1*10^4; %test specific conductance for membrane S/cm^2
r = .5; %test radius for neurite in micrometers
Ri = 100; %test axial resistance Ohm * cm

rCm = r / 10^4; %convert micrometers to centimeters
Rm = 1 / Gm; %Calculate specific membrane resistance from specific membrane conductance
rm = Rm / (2 * pi * rCm); % Calculate membrane resistance from specific membrane resitance and radius
ra = Ri / (pi * rCm^2); % Calculate axial resistance from specific axial resistance and radius
lcCM = sqrt(rm/ra); % Calculate length constant (cm) from membrane resistance and axial resistance
lc = lcCM * 10^4 % Convert centimeters to micrometers











