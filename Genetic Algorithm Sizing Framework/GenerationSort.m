
numDV = 2;


% load Generation.mat
load('Generation.mat')

% 'c' is the number of generations evaluated
NumGen = c;



% start a for loop running over each generation
for i = 1:NumGen
   

    % form generation name
    GenName = sprintf('Generation%0.0f',i-1);

    % load generation
    load(GenName)

    % create a folder for this generation
    mkdir(GenName)

    % variable 'b' contains the design variable and objective function info for
    % the population of this generation
    DV = b(:,1:numDV);

    

    % run over each population member
    for j = 1:1:height(b)

        % form vehicle file name
        VehName0 =  strrep(sprintf('Aircraft_WL%0.4f_AR%0.4f',DV(j,1),DV(j,2)),'.','_');
        VehName = sprintf('%s.mat',VehName0);
        try
        % move the file into the correct folder
        copyfile(VehName,sprintf('%s/',GenName));
        catch ERRMSG
        end

    end

    % create a copy of the generation file into this folder 
    copyfile(sprintf('%s.mat',GenName),sprintf('%s/',GenName))



end