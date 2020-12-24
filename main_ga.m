clear all
clc
close all

addpath('time');
addpath('functions');
%% Data

%Planetary constants [km^3/s^2]
mu_nep   = astroConstants(18);      %Neptune  
mu_earth = astroConstants(13);      %Earth
mu_merc  = astroConstants(11);      %Mercury



%Radius and height of atmosphere of Earth
Re   = astroConstants(23);        %[km]
hatm = 100;                       %[km] 

%Earliest departure and latest arrival dates
dep = date2mjd2000([2027,04,1,0,0,0]);
arr = date2mjd2000([2067,04,1,0,0,0]);

%Orbital parameters at earliest departure date
%[a,e,i,w,W,f]
% a [km]    Semi-major axis
% e [-]     Eccentricity
% i [rad]   Inclination
% w [rad]   Argument of perigee
% W [rad]   Right Ascension of Ascending Node
% f [rad]   True anomaly

[kep_nep,ksun] = uplanet(dep, 8);       %Neptune
[kep_earth,ksun] = uplanet(dep,3 );     %Earth
[kep_merc,ksun] = uplanet(dep, 1);      %Mercury

% Synodic periods

[Tsyn_ne,Tsyn_me,Tsyn_nm,ris] = synodic_periods(kep_nep,kep_earth,kep_merc,ksun);




%% Hohmann tranfers tof

[tof_h_nep_earth,tof_h_earth_merc] = tofs_hohmann(dep,kep_nep,kep_earth,kep_merc,ksun);

%Minimum and maximum departure and arrival times
dep2h = arr - tof_h_nep_earth - tof_h_earth_merc;       %Latest departure date based on results of Hohmann transfer
arrh_min_earth = dep + tof_h_nep_earth;                 %Earliest arrival date on planet 2 based on results of Hohmann transfer
arrh_max_earth = dep2h + tof_h_nep_earth;               %Latest arrival date on planet 2 based on results of Hohmann transfer
arrh_min = dep + tof_h_nep_earth + tof_h_earth_merc;    %Earliest arrival date on planet 3 based on results of Hohmann transfer


%% Time windows of time for the 3 phases
n = 15;

departure_vec = linspace(dep,dep2h,n); 
arrival_vec_earth = linspace(arrh_min_earth,arrh_max_earth,n);
arrival_vec = linspace(arrh_min,arr,n);

% departure_vec = [dep:200:arr];                            
% arrival_vec_earth = [dep:200:arr];
% arrival_vec = [dep:200:arr];

A = [-1,0,0;1,0,0;1,1,0;1,1,1];
b = [-dep;dep2h;arrh_max_earth;arr];
[X1,FVAL1] = ga(@sum_ga,3,A,b);

lungh = (dep2h-dep)/20; %days

vec = [dep:lungh:dep2h];

for i = 1:length(vec)-1
    A = [-1,0,0; 1,0,0; 1,1,0; 1,1,1];
    b = [-vec(i); vec(i+1); arrh_max_earth; arr];
    [X(i,:),FVAL(i)] = ga(@sum_ga,3,A,b);
end























