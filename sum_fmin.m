function deltav_fmin= sum_fmin(x)

% x vettore che ha dep,TOF1,TOF2

Re   = astroConstants(23);        %[km]
hatm = 100;                       %[km] 
mu_earth = astroConstants(13);      %Earth



if x(2) > 0 && x(3) > 0 && x(3) < 2000 

    %Departure ephemeris
    [kep_nep,ksun] = uplanet(x(1), 8);
    [r1,v1] =kep2car(kep_nep,ksun);   
   

    %Arrival ephemeris
    [kep_merc,ksun] = uplanet(x(1)+x(2)+x(3), 1);
    [r2,v2] =kep2car(kep_merc,ksun);   


    %Fly_by ephemeris
    [kep_earth,ksun] = uplanet(x(1)+x(2), 3);
    [r_m,v_m] =kep2car(kep_earth,ksun);
            
    Vpl = v_m;         %Velocity of planet 2 (Earth) in heliocentric frame
            
    %Dv for heliocentric leg from Neptune to Earth
    [~,delta_v1,~,VI,VF] = comp_dv_fmin(r1,v1,r_m,v_m,x(1),x(1)+x(2),ksun);            
    V_M = VF';            
            
    %Dv for heliocentric leg from Earth to Mercury
    [~,~,delta_v3,VI,VF] = comp_dv_fmin(r_m,v_m,r2,v2,x(1)+x(2),x(1)+x(2)+x(3),ksun);
    V_P = VI';        
            
    %Entry and exit velocities in SOI of Earth
    vinfM = V_M - Vpl;        
    vinfP = V_P - Vpl; 
            
            
    [deltav_perig,rp] = flybyPow_fmin(vinfM,vinfP,mu_earth,hatm,Re+hatm); 
    if rp <= Re+hatm
         deltav_perig = 40;
    end
    deltav_fmin = delta_v1+delta_v3+deltav_perig;

else
    deltav_fmin = 40;
end
                
          
            
                

             


