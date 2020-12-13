function [deltat_hom_days,deltat_hom_days2] = tofs_hohmann(dep_min,kep_nep,kep_earth,kep_merc,ksun)



%time of flight of Homhann tf from nep to earth
a_hom = (kep_nep(1)+kep_earth(1))/2;
deltat_hom = pi*sqrt(a_hom^3/ksun);
deltat_hom_days = deltat_hom/(3600*24);

%time of flight of Homhann tf from earth to mercury
a_hom2 = (kep_earth(1)+kep_merc(1))/2;
deltat_hom2 = pi*sqrt(a_hom2^3/ksun);
deltat_hom_days2 = deltat_hom2/(3600*24);

% dep2 = arr_max - deltat_hom_days - deltat_hom_days2;
% arr_min_earth = dep_min + deltat_hom_days;
% arr_max_earth = dep2 + deltat_hom_days;