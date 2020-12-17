function [Tsyn_ne,Tsyn_me,Tsyn_nm,ris] = synodic_periods(kep_nep,kep_earth,kep_merc,ksun)

%Orbital period [s]
T_nep = 2*pi*sqrt(kep_nep(1)^3/ksun);
T_earth = 2*pi*sqrt(kep_earth(1)^3/ksun);
T_merc = 2*pi*sqrt(kep_merc(1)^3/ksun);

%Synodic periods [s]
Tsyn_ne = T_nep*T_earth/abs(T_nep-T_earth);
Tsyn_me = T_merc*T_earth/abs(T_merc-T_earth);
Tsyn_nm = T_nep*T_merc/abs(T_nep-T_merc);

l = Tsyn_ne/Tsyn_me;

a = Tsyn_ne*35;
b = Tsyn_me*111;
ris=a-b;
