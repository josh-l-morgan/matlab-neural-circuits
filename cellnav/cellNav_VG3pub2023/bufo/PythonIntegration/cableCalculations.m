
%%When x = L (distance = length constant), V(x) = 0.63 Vm (maximum voltage)
%%mammalian nerve fiber of 1 micron diam has Length constant of about 0.2
%%mm



Sm = 0.001; % specific membrane conductance
Ra = 100 % internal resistivity, Ohms * cm
ro = 0 % external resistance
rad = .2 % radius im um
Cm = 1; % specific membrane capacitance = uF/cm^2

Rm = 1/Sm;% %membrane resistivity, Ohms * cm^2
r = rad/10^5; % radius in centimeters
rm = Rm/ (2 * pi * r) % membrane resitance, Ohms * cm
ri = Ra / (pi * r^2) % internal resistance, Ohms / cm


L = sqrt(rm/ri)
L2 = sqrt(r * Rm / (2 * Ra))
tau = Rm * Cm

Vm = 1; %milivolt
I = Vm/Rm


Lu = L * 10^5 %length constant in um
x = [0:(Lu*100)]/10;
Vrise = Vm * exp(-x/Lu);
Vfall = Vm * (1-exp(-x/Lu));

plot(x,Vrise,'b')
hold on
plot(x,Vfall,'r')
hold off
Vrise(round(Lu*10)+1)


%% From Koch Segeve Methods in Neuronal modeling
lambda=sqrt(Rm/Ra*r/2) * 10^5 %length constant in um




