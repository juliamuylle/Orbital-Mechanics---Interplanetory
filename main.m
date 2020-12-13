% dep = nep
% fb = earth
% ar = Mercr
% early 1/04/2027
% lat 1/04/2067

clear all
clc
close all
addpath('time');

mu_nep = astroConstants(18);
mu_earth = astroConstants(13);
mu_merc = astroConstants(11);
Re = astroConstants(23);        %[km] Radius of Earth

dep = date2mjd2000([2027,04,1,0,0,0]);
arr = date2mjd2000([2067,04,1,0,0,0]);
[kep_nep,ksun] = uplanet(dep, 8);
[kep_earth,ksun] = uplanet(dep,3 );
[kep_merc,ksun] = uplanet(dep, 1);

T_nep = 2*pi*sqrt(kep_nep(1)^3/ksun);
T_earth = 2*pi*sqrt(kep_earth(1)^3/ksun);
T_merc = 2*pi*sqrt(kep_merc(1)^3/ksun);

Tsyn_ne = T_nep*T_earth/abs(T_nep-T_earth);
Tsyn_me = T_merc*T_earth/abs(T_merc-T_earth);
Tsyn_nm = T_nep*T_merc/abs(T_nep-T_merc);

j = Tsyn_ne/Tsyn_me;

a = Tsyn_ne*11;
b = Tsyn_me*35;
ris = a-b;

% a = Tsyn_ne*35;
% b = Tsyn_me*111;
% a-b



%% Hohmann tranfers

[tof_h_nep_earth,tof_h_eart_merc] = tofs_hohmann(dep,kep_nep,kep_earth,kep_merc,ksun);



%% transfer from neptune to earth
dep2h = arr - tof_h_nep_earth - tof_h_eart_merc;
arrh_min_earth = dep + tof_h_nep_earth;
arrh_max_earth = dep2h + tof_h_nep_earth;
arrh_min = dep + tof_h_nep_earth + tof_h_eart_merc;

departure_vec = [dep:100:dep2h];
arrival_vec_earth = [arrh_min_earth:100:arrh_max_earth];
arrival_vec = [arrh_min:100:arr];



%% fly-by

counti = 1;
for i = departure_vec
   [kep_nep,ksun] = uplanet(i, 8);
   [r1(:,counti),v1(:,counti)] =kep2car(kep_nep,ksun);

   countj = 1;
   for j =arrival_vec
        [kep_merc,ksun] = uplanet(j, 1);
        [r2(:,countj),v2(:,countj)] =kep2car(kep_merc,ksun); 
         
         countk = 1;
         for k = arrival_vec_earth
            [kep_earth,ksun] = uplanet(k, 3);
            [r_m(:,countk),v_m(:,countk)] =kep2car(kep_earth,ksun);
            
            
            
            Vpl = v_m(:,countk);
            
            %dv for nep_earth
            [delta_v,VI,VF] = comp_dv(r1,v1,r_m,v_m,departure_vec,arrival_vec_earth,ksun);            
            V_M = VF;
            
            %dv for earth_merc
            [delta_v,VI,VF] = comp_dv(r_m,v_m,r2,v2,departure_vec,arrival_vec_earth,ksun);
            V_P = VI;
            
            vinfM = V_M-Vpl;        
            vinfP = V_P-Vpl; 
        
            [deltav_perig,rp,delta,arcs] = flybyPow(vinfM,vinfP,mu_earth,Re);
 
            %Total cost
            deltaV_tot = V_P-V_M;
            
            
            countk = countk+1;
         end
    
     countj = countj+1;
   end
    counti = counti+1;
end

%dv for nep_earth
delta_v1 = comp_dv(r1,v1,r_m,v_m,departure_vec,arrival_vec_earth,ksun);


%dv for earth_merc
delta_v3 = comp_dv(r_m,v_m,r2,v2,departure_vec,arrival_vec_earth,ksun);


delta_v = delta_v1+delta_v3;
deltavmin = min(delta_v);
deltavmin = min(deltavmin);

%%

for k1 = 1:length(departure_vec)
    t1_plot(k1) = datenum(mjd20002date(departure_vec(k1)));
end
for k2 = 1:length(arrival_vec)
    t2_plot(k2) = datenum(mjd20002date(arrival_vec(k2)));
end

figure(1)
[c,h]=contour(t1_plot,t2_plot,delta_v',floor(deltavmin)+(0:1:10),'Showtext','off');
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



