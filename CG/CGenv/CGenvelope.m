% Calculating CG envelope
% Prof. Chakraborty - 
% fall 2025
%--------------------------------------------------------------------------
clc; clear; close all

locweight =  1.0e+03 * [0.0055    0.5000; ...
                            0.0090    0.3000;
                            0.0130    2.0000;
                            0.0200    9.8062;
                            0.0300    9.8062;
                            0.0436    2.4500;
                            0.0225    1.3047;
                            0.0278    5.4097;
                            0.0428    0.3968;
                            0.0408    0.5490;
                            0.0268    2.0300;
                            0.0322    2.0300;
                            0.0141    0.3760];

beta = [0.9985 0.9756 0.9930 0.8932 0.9525 0.9710 0.9930 0.8643 0.9443 0.9764 0.9974]';
    
cginit = (sum(locweight(:,1).*locweight(:,2)))./sum(locweight(:,2))

for i = 1:length(beta)
    
    weight = locweight;

    fburn = sum(locweight(:,2))*(1-prod(beta(1:i)));
    
    if (locweight(5,2)-fburn>0)

        weight(5,2) = weight(5,2)-fburn;
        CG(i)=(sum(weight(:,1).*weight(:,2)))./sum(weight(:,2));

    else

        weight(4,2) = weight(4,2) - fburn + weight(5,2);
        weight(5,2) = 0;
        CG(i)=(sum(weight(:,1).*weight(:,2)))./sum(weight(:,2));

    end


    

end

CG = [26.3203 CG];
miss = 12:-1:1;
plot(CG,miss,"-o")
grid on
xlabel("CG location (ft)")

fprintf("inital weight: " + sum(locweight(4:5,2))+" (lbs).\n")

fprintf("leftover fuel: " + string(sum(weight(4:5,2)))+" (lbs).\n\n")