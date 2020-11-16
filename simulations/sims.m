nWords_all = [20 100 1000 10000];
actioncost_all = [2 5 10 25 50 100];

for nWords_ind = 1
    for actioncost_ind = 1:length(actioncost_all)
        % the four parameters we manipulate
        nWords = nWords_all(nWords_ind);
        actioncost = actioncost_all(actioncost_ind);
        %nToEval = [1:20 40 80 nWords];
        nToEval = [1:20];
        r = [0 .25 .5 .75 1];
        
        nEnvs = length(r);
        variance = 1;
        mu = [0 0];
        nSizes = length(nToEval);
        
        nAgents = 10000;
        mfearnings = zeros(nEnvs, nSizes, nAgents);
        mbearnings = zeros(nEnvs, nSizes, nAgents);
        csearnings = zeros(nEnvs, nSizes, nAgents);
        rcsearnings = zeros(nEnvs, nSizes, nAgents);
        
        parfor env = 1:nEnvs
            cur_r = r(env);
            
            sigma = [variance cur_r * variance; cur_r * variance variance];
            
            for curSize = 1:nSizes
                cur_nToEval = nToEval(curSize);
                
                for agent = 1:nAgents
                    temp = mvnrnd(mu, sigma, nWords);
                    s1re = temp(:, 1);
                    s2re = temp(:, 2);
                    
                    % run cs
                    [~, options] = sort(s1re, 'descend');
                    toEval = options(1:cur_nToEval);
                    [~, choice_ind] = max(s2re(toEval));
                    choice = toEval(choice_ind);
                    csearnings(env, curSize, agent) = s2re(choice);
                end
            end
        end
        
        %% Do optimality computation, save data to be graphed
        
        cs_mean = mean(csearnings, 3);
        output = cs_mean';
        
        % get # of decisions possible for each cs size
        budget = 1000 + actioncost;
        
        nDecisions = [budget / actioncost, budget ./ (actioncost + nToEval(2:end))];
        
        % Compute the average decision value
        decisionVals = bsxfun(@times, nDecisions', output);
        decisionVals_nonzero = decisionVals(:, 2:end); % value-guided CS sampling
        decisionVals_zero = repmat(decisionVals(:, 1), 1, size(decisionVals_nonzero, 2)); % random CS sampling
        maxVals = max(decisionVals_nonzero);
        
        % normalize
        normalizedVals_nonzero = bsxfun(@rdivide, decisionVals_nonzero, maxVals);
        normalizedVals_zero = bsxfun(@rdivide, decisionVals_zero, maxVals);
        
        normVals = [normalizedVals_nonzero normalizedVals_zero];
        output2 = [normVals(:), repmat(repelem(r(2:end)', nSizes), 2, 1), repmat(nToEval', (nEnvs - 1) * 2, 1), repelem([0 1], nSizes * (nEnvs - 1))'];
        
        csvwrite(['data/' num2str(nWords) '_' num2str(actioncost) '.csv'], output2);
    end
end