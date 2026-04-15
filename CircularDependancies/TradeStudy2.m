clearvars;
clc;

displayResults = false;

valuesPerInput = 2;

WL = linspace(80, 105, valuesPerInput);
AR = linspace(2.5, 5, valuesPerInput);
TR = linspace(0.2, 0.55, valuesPerInput);
SA = linspace(25, 55, valuesPerInput);

results  = cell(valuesPerInput,valuesPerInput,valuesPerInput,valuesPerInput);
MTOWs = cell(valuesPerInput,valuesPerInput,valuesPerInput,valuesPerInput);

for i = 1:valuesPerInput
    for j = 1:valuesPerInput
        for k = 1:valuesPerInput
            for l = 1:valuesPerInput
                clear aircraft;

                aircraft = struct;
                aircraft.constants.wingLoading = WL(i);
                aircraft.wing.AR = AR(j);
                aircraft.wing.taper_ratio = TR(k);
                aircraft.wing.sweep = SA(l);

                try 
                    aircraft = Main_forTS(aircraft);
                    results{i,j,k,l} = aircraft;
                    MTOWs{i,j,k,l} = aircraft.weight.total;
                catch ERRMSG
                    results{i,j,k,l} = ERRMSG.message;
                    fprintf("Failure: ")
                end

                fprintf("i: %u, j: %u, k: %u, l: %u\n",i,j,k,l);


            end % l
        end % k
    end % j 
end % i