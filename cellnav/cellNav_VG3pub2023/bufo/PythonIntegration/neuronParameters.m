function[p] = neuronParameters(pick)
%%Default Neuron parameters

if pick == 0
    p = {'Long Length Constant','Default', 'AII 1' 'AMC 1'};

elseif pick == 1
    %%Neuron default
    p.cm = 1;
    p.Ra = 1;
    p.Rm = 0.001;
    
elseif pic == 2
        %%Neuron default
    p.cm = 1;
    p.Ra = 100;
    p.Rm = 0.001;

elseif pick == 3
    %%Electrotonic signal processing in AII amacrine cells: compartmentalmodels
    %%and passive membrane properties for a gap junction-coupled retinal neuron
    %%Bas‑Jan Zandt Margaret Lin Veruki Espen Hartveit
    %%Brain Structure and Function (2018) 223:3383–3410

    p.Ra = 200; % (Ri) Oms cm, cytoplasmic resistivity
    p.cm = 0.91; %uF cm^-2, specific membran capacitance
    p.Rm = .00302; %Oms cm^2, specific membrane resitantce
    %p.Rs = 25; % MOhm, series resistance,

elseif pick == 4
    %%Inward rectifying currents stabilize the membrane potential in
    %%dendrites of mouse amacrine cells: Patch-clamp recordings and
    %%single-cell RT-PCR
    %%Amane Koizumi, Tatjana C. Jakobs, Richard H. Masland
   
    p.cm = 1;
    p.Ra = 150; %Ohms cm
    membraneConductance = 3*10^5; %g/cm^2
    p.Rm = 1/membraneConductance ;%?

else
    
    p.cm = 1;
    p.Ra = 100;
    p.Rm = 0.001;

end