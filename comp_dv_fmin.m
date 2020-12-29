function [delta_v,delta_v1,delta_v2,VI,VF] = comp_dv_fmin(r1,v1,r2,v2,dep,ar,ksun)


orbitType = 0; %0 if prog; 1 if retrog
Nrev = 0; 
Ncase = 0;
optionsLMR = 0;

        TOF = (ar-dep)*24*3600;
        
            [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1,r2,TOF,ksun,orbitType,Nrev,Ncase,optionsLMR);
            delta_v1 = norm(VI'- v1);
            delta_v2 = norm(v2-VF');
            delta_v = delta_v1+delta_v2;
        