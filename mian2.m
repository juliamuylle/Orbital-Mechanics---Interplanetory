% dep = nep
% fb = earth
% ar = Mercr
% earliest 1/04/2027
% lat 1/04/2067

clear all
clc
close all

addpath('time');
%% Data

%Planetary constants [km^3/s^2]
mu_nep = astroConstants(18);        %Neptune  
mu_earth = astroConstants(13);      %Earth
mu_merc = astroConstants(11);       %Mercury

hatm = 100;

%Radius of Earth
Re = astroConstants(23);        %[km]

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

%Orbital period [s]
T_nep = 2*pi*sqrt(kep_nep(1)^3/ksun);
T_earth = 2*pi*sqrt(kep_earth(1)^3/ksun);
T_merc = 2*pi*sqrt(kep_merc(1)^3/ksun);

%Synodic periods [s]
Tsyn_ne = T_nep*T_earth/abs(T_nep-T_earth);
Tsyn_me = T_merc*T_earth/abs(T_merc-T_earth);
Tsyn_nm = T_nep*T_merc/abs(T_nep-T_merc);

l = Tsyn_ne/Tsyn_me;

% a = Tsyn_ne*11;
% b = Tsyn_me*35;
% ris = a-b;

a = Tsyn_ne*35;
b = Tsyn_me*111;
ris=a-b;


%% Hohmann tranfers tof

[tof_h_nep_earth,tof_h_earth_merc] = tofs_hohmann(dep,kep_nep,kep_earth,kep_merc,ksun);


%% Transfer from Neptune to Earth
% dep2h = arr - tof_h_nep_earth - tof_h_earth_merc;       %Latest departure date based on results of Hohmann transfer
% arrh_min_earth = dep + tof_h_nep_earth;                 %Earliest arrival date on planet 2 based on results of Hohmann transfer
% arrh_max_earth = dep2h + tof_h_nep_earth;               %Latest arrival date on planet 2 based on results of Hohmann transfer
% arrh_min = dep + tof_h_nep_earth + tof_h_earth_merc;    %Earliest arrival date on planet 3 based on results of Hohmann transfer

%Windows of time for the 3 phases
% departure_vec = [dep:500:dep2h];                            
% arrival_vec_earth = [arrh_min_earth:500:arrh_max_earth];
% arrival_vec = [arrh_min:500:arr];

departure_vec = [dep:500:arr];                            
arrival_vec_earth = [dep:500:arr];
arrival_vec = [dep:500:arr];

%% 

DV1=[];
DV3=[];

x = [Re+hatm:1:50000];

%Departure ephemeris

for i = 1 : length(departure_vec)
   [kep_nep,ksun] = uplanet(departure_vec(i), 8);
   [r1,v1] =kep2car(kep_nep,ksun);
   r1_vec(:,i) = r1;
   v1_vec(:,i) = v1;

    
%Arrival ephemeris

    
   for j = 1 : length(arrival_vec)
        [kep_merc,ksun] = uplanet(arrival_vec(j), 1);
        [r2,v2] =kep2car(kep_merc,ksun); 
         r2_vec(:,i) = r2;
         v2_vec(:,i) = v2;


%Fly_by ephemeris         

    
         for k = 1 : length(arrival_vec_earth)
            [kep_earth,ksun] = uplanet(arrival_vec_earth(k), 3);
            [r_m,v_m] =kep2car(kep_earth,ksun);
            Vpl = v_m;         %Velocity of planet 2 (Earth) in heliocentric frame
            
            %Dv for heliocentric leg from Neptune to Earth
            [delta_v1,VI,VF] = comp_dv(r1,v1,r_m,v_m,departure_vec(i),arrival_vec_earth(k),ksun);            
            V_M = VF';
            DV1=[DV1;delta_v1];
            
            %Dv for heliocentric leg from Earth to Mercury
            [delta_v3,VI,VF] = comp_dv(r_m,v_m,r2,v2,arrival_vec_earth(k),arrival_vec(j),ksun);
            V_P = VI';
            DV3=[DV3;delta_v3];
            
            %Entry and exit velocities in SOI of Earth
            vinfM = V_M-Vpl;        
            vinfP = V_P-Vpl; 
            
            if arrival_vec_earth(k)-departure_vec(i) > 0 && arrival_vec(j) - arrival_vec_earth(k) > 0
                [deltav_perig(i,j,k),rp(i,j,k),delta,arcs,FUN,delta_true,deltaM,deltaP,eM,eP] = flybyPow(vinfM,vinfP,mu_earth,hatm,Re+hatm); 

            else
                rp(i,j,k) = 0;
                deltav_perig(i,j,k) = NaN;
            end
                %             
%epp = eP(x);
%             emm = eM(x);
%             figure(15)
%             plot(x,eP(x))
%             hold on
%   
%         
%             figure(10)
%             plot(x,deltaM(x)/3.14*180)
%             plot(x,delta_true*ones(length(x),1)/3.14*180)
%             hold on

             
            
         end
   end
end

%[deltav_perig,rp,delta,arcs,FUN,delta_true,eM] = flybyPow(vinfM,vinfP,mu_earth,hatm,Re+hatm); 
% FUN(307)
% rp = [Re+hatm:1:50000];
% 
% figure(1)
% plot(rp,eM(rp))

%Total cost  ?????
deltaV_tot = delta_v1+delta_v3;


deltavmin = min(deltaV_tot);
deltavmin = min(deltavmin);

%%
%Dates in Matlab format
for k1 = 1:length(departure_vec)
    t1_plot(k1) = datenum(mjd20002date(departure_vec(k1)));
end
for k2 = 1:length(arrival_vec)
    t2_plot(k2) = datenum(mjd20002date(arrival_vec(k2)));
end

figure(1)
[c,h]=contour(t1_plot,t2_plot,deltaV_tot',floor(deltavmin)+(0:1:10),'Showtext','off');
%clabel(C,h,floor(deltavmin)+[0:1:5])
caxis(floor(deltavmin)+[0 10]);
caxis('manual')
xtickangle(45) % put label in a inclined angle
ytickangle(45)
datetick('x','yyyy mm dd','keeplimits')
datetick('y','yyyy mm dd','keeplimits')
ylabel('Arrival date')
xlabel('Departure date')
hcb = colorbar
hcb.Title.String = '$Delta v$ [km/s]';
hcb.Title.Interpreter = 'latex';
%  ha = gca;
% ha.FontSize = 13;
% hcb.Title.FontSize = 15;

%Constant tof lines
hold on
[c2,h2] = contour(t1_plot,t2_plot, t2_plot'-t1_plot,[60,120,180,240,300],'k');
clabel(c2,h2); %give the name to line

% DVtot(DVtot>6) = NaN %I take only the elem of DVtot that are larger than
% 6 and I out them = NaN
%for es 4 for ex I can do DVtot(DV1<vinf) = NaN to eliminate the elem of
%DVtot in correspondance of the ones of DV1 that are higher than the
%constraint NB the 2 matrices must have same size

%find delta vmin, date of arrival and departure

[row,col] = find(delta_v==deltavmin);
departure = departure_vec(row);
arrival = arrival_vec(col);


%propagate the orbit for the mission with dv min
r1_min = r1(:,row);
v1_min = v1(:,row);
r2_min = r2(:,col);
v2_min = v2(:,col);
MU = ksun;
TOF = (arrival-departure)*24*3600;
orbitType = 0; %0 if prog; 1 if retrog
Nrev = 0; 
Ncase = 0;
optionsLMR = 0;
[A,P,E,ERROR,VI,VF,~,~] = lambertMR(r1_min,r2_min,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);

T = 2*pi*sqrt(A^3/MU);
tspan = [0 T];
options = odeset ( 'RelTol', 1e-13,'AbsTol', 1e-14 );


y01 = [r1_min,v1_min];
[t1,Y1] = ode113(@(t,y) odefun(MU,y,t),tspan,y01,options);

y0T = [r1_min',VI];
[tT,YT] = ode113(@(t,y) odefun(MU,y,t),tspan,y0T,options);

y02 = [r2_min,v2_min];
[t2,Y2] = ode113(@(t,y) odefun(MU,y,t),tspan,y02,options);


figure(2)
plot (Y1(:,1),Y1(:,2))
hold on
plot(YT(:,1),YT(:,2))
hold on
plot (Y2(:,1),Y2(:,2))
hold on
plot(r1_min(1),r1_min(2),'*')
hold on
plot(r2_min(1),r2_min(2),'*')
grid on
xlabel('rx')
ylabel('ry')
legend('Orbit1','Tranfer arc','Orbit2','P1','P2')
title('Orbits of the manoeuver')
axis equal
