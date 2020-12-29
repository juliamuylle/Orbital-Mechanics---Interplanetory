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
arr = date2mjd2000([2064,08,1,0,0,0]);


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
dep2h = date2mjd2000([2033,01,1,0,0,0]);       %Latest departure date based on results of Hohmann transfer
arrh_min_earth = date2mjd2000([2058,01,1,0,0,0]);                 %Earliest arrival date on planet 2 based on results of Hohmann transfer
arrh_max_earth =date2mjd2000([2061,01,1,0,0,0]);              %Latest arrival date on planet 2 based on results of Hohmann transfer
arrh_min = date2mjd2000([2061,01,2,0,0,0]); 

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


lungh = (dep2h-dep)/300; %days

vec = [dep:lungh:dep2h];



for i = 1:length(vec)-1
    A = [-1,0,0; 1,0,0; 1,1,0; 1,1,1];
    b = [-vec(i); vec(i+1); arrh_max_earth; arr];   
    x0 = [vec(i)+(vec(i)+vec(i+1))/2; arrh_min_earth+(arrh_max_earth-arrh_min_earth)/2-(vec(i)-vec(i+1))/2-vec(i); arrh_min+(arr-arrh_min)/2-(arrh_max_earth-arrh_min_earth)/2-arrh_min_earth];
    [X(i,:),FVAL(i)] = fmincon(@sum_fmin,x0,A,b);
end























