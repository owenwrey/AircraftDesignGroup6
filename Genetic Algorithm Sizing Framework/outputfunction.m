function [state, options, optchanged] = outputfunction(options, state, ~)
    optchanged = false;

    a = state.Generation;
    b = [state.Population state.Score];

    % Save generation data
    SaveFileName = ['Generation', num2str(a), '.mat'];
    save(SaveFileName, 'b')

    % Stop if population has fully collapsed
    spread = max(state.Score) - min(state.Score);
    if max(abs(spread)) < 1e-7
        state.StopFlag = 'Population converged';
    end
end