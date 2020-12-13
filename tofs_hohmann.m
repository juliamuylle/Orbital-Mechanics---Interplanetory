function [deltat_hom_days,deltat_hom_days2] = tofs_hohmann(dep_min,kep_pl1,kep_pl2,kep_pl3,ksun)

%tofs_hohmann returns the time of flight expressed in days of Hohmann transfer from planet 1
%to planet 2 and from planet 2 to planet 3
%
%PROTOTYPE: 
%     [deltat_hom_days,deltat_hom_days2] = tofs_hohmann(dep_min,kep_nep,kep_earth,kep_merc,ksun)
% 
% INPUT:
%     dep_min [1]           Earliest departure date
%     kep_pl1 [1x6]     	Keplerian parameters of planet 1 at dep_min
%     kep_pl2 [1x6]     	Keplerian parameters of planet 2 at dep_min
%     kep_pl3 [1x6]         Keplerian parameters of planet 3 at dep_min
%     ksun [1]              Sun planetary constant [km^3/s^2]
%     
% OUTPUT:
%     deltat_hom_days [1]   Time of flight of Homhann transfer from Neptune
%                           to Earth [days]
%     deltat_hom_days2 [1]  Time of flight of Homhann transfer from Earth
%                           to Mercury [days]
%
% CONTRIBUTORS
%       Bertolini Edoardo
%       Busi Silvia
%       Muylle Julia
%       Pellegrini Matias
%
% VERSIONS
%
% 13/12/2020: First Version

%Time of flight of Homhann transfer from Neptune to Earth
a_hom = (kep_pl1(1)+kep_pl2(1))/2;          %Semi-major axis [km]
deltat_hom = pi*sqrt(a_hom^3/ksun);         %Tof [s]
deltat_hom_days = deltat_hom/(3600*24);     %Tof[days]

%Time of flight of Homhann transfer from Earth to Mercury
a_hom2 = (kep_pl2(1)+kep_pl3(1))/2;         %Semi-major axis [km]
deltat_hom2 = pi*sqrt(a_hom2^3/ksun);       %Tof [s]
deltat_hom_days2 = deltat_hom2/(3600*24);   %Tof [days]

% dep2 = arr_max - deltat_hom_days - deltat_hom_days2;
% arr_min_earth = dep_min + deltat_hom_days;
% arr_max_earth = dep2 + deltat_hom_days;
