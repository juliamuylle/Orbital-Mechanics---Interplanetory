function [delta_v,VI,VF] = comp_dv(r1,v1,r2,v2,dep,ar,ksun)




orbitType = 0; %0 if prog; 1 if retrog
Nrev = 0; 
Ncase = 0;
optionsLMR = 0;

for i = 1:length(r1(1,:))
    for j = 1:length(r2(1,:))
        TOF = (ar(j)-dep(i))*24*3600;
        if TOF > 0
            [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1(:,i),r2(:,j),TOF,ksun,orbitType,Nrev,Ncase,optionsLMR);
            delta_v1 =VI'- v1(:,i);
            delta_v2 = v2(:,j)-VF';
            delta_v(i,j) = norm(delta_v1)+norm(delta_v2);
        else
            delta_v(i,j) = NaN;
        end
    end
end
        
    

