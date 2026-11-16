%% ========================================================================
%  pcscore_leave3_dmel.m
%
%  INPUT-RESILIENCE ANALYSIS (leave-THREE-fruits-out PC score) — D. melanogaster
%
%  Purpose
%  -------
%  Same idea as the leave-one-fruit-out analysis for noni/pandan, but for the
%  larger D. melanogaster dataset (2,491 odorants x 33 fruits). Because leaving
%  3 fruits out gives C(33,3) = 5,456 combinations, we do not enumerate them all
%  by default: we draw a fixed number of RANDOM 3-fruit combinations (with a set
%  random seed for reproducibility) and recompute the PC-score pipeline for each.
%  The spread of PC scores across these iterations shows how robust each
%  odorant's ranking is to which fruits are sampled.
%
%  The PC-score computation is identical to the main analyses
%  (chem_pca_plus_PCscore.m): within-fruit (row-wise) z-scoring, PCA with fruits
%  as observations and odorants as variables, PC score = sum over the first three
%  PCs of |loading| * percent variance explained.
%
%  Input
%  -----
%  CSV with odorants as ROWS, fruits as COLUMNS:
%     - column 1  : odorant names
%     - columns 2+: one column per fruit
%
%  Output
%  ------
%  One figure: x-axis = the TOP-N odorants by median PC score (the full 2,491
%  cannot be shown legibly; the analysis still uses all odorants, only the
%  display is truncated to the top ranked). y-axis = PC score, one grey dot per
%  leave-3-out iteration, over a grey boxplot.
%  ========================================================================

clear; close all; clc;

%% ------------------------- USER SETTINGS --------------------------------
dmelFile   = 'dmel_data.csv';

nPCs       = 3;      % leading PCs used in the PC score
leaveOut   = 3;      % fruits removed per iteration
nIterWanted = 1000;  % number of random leave-3-out combinations to draw
                     %   (if this exceeds the number of unique combinations,
                     %    the script falls back to exhaustive enumeration)
topN       = 10;     % number of top odorants (by median PC score) to display
randSeed   = 1;      % fixed seed for reproducible figure

%% ------------------------- RUN ------------------------------------------
runResilienceRandom(dmelFile, 'D. melanogaster', nPCs, leaveOut, ...
                    nIterWanted, topN, randSeed);


%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function runResilienceRandom(csvFile, datasetName, nPCs, leaveOut, ...
                             nIterWanted, topN, randSeed)

    %% ---- Load data -----------------------------------------------------
    T = readtable(csvFile, 'VariableNamingRule','preserve');
    odourNames = strtrim(string(T{:,1}));
    X          = T{:,2:end};                 % odours x fruits
    if ~isnumeric(X)
        X = str2double(string(X));
    end

    nOdours = size(X,1);
    nFruits = size(X,2);
    fprintf('%s: %d odorants x %d fruits\n', datasetName, nOdours, nFruits);

    %% ---- Build the set of leave-out combinations -----------------------
    % Each row of "combos" lists the fruit indices REMOVED in that iteration.
    rng(randSeed);                               % reproducible sampling
    maxCombos = nchoosek(nFruits, leaveOut);

    if nIterWanted >= maxCombos
        % Few enough to enumerate every combination exactly
        combos = nchoosek(1:nFruits, leaveOut);
        fprintf(['Enumerating all %d leave-%d-out combinations ' ...
                 '(<= requested %d).\n'], maxCombos, leaveOut, nIterWanted);
    else
        % Draw unique random combinations without replacement
        combos = zeros(nIterWanted, leaveOut);
        seen   = containers.Map('KeyType','char','ValueType','logical');
        count  = 0;
        while count < nIterWanted
            c   = sort(randperm(nFruits, leaveOut));   % one random combo
            key = sprintf('%d-', c);
            if ~isKey(seen, key)
                seen(key)      = true;
                count          = count + 1;
                combos(count,:) = c;
            end
        end
        fprintf('Drew %d unique random leave-%d-out combinations (of %d possible).\n', ...
                nIterWanted, leaveOut, maxCombos);
    end
    nIter = size(combos,1);

    %% ---- Compute PC score for each iteration ---------------------------
    PCscoreAll = nan(nOdours, nIter);
    for it = 1:nIter
        keep = true(1, nFruits);
        keep(combos(it,:)) = false;
        PCscoreAll(:,it) = computePCscore(X(:,keep), nPCs);
    end

    %% ---- Rank odours by median PC score, keep top-N for display --------
    medScore   = median(PCscoreAll, 2, 'omitnan');
    [~, order] = sort(medScore, 'descend');

    keepIdx    = order(1:min(topN, nOdours));
    PCscoreTop = PCscoreAll(keepIdx, :);
    namesTop   = odourNames(keepIdx);

    %% ---- Plot ----------------------------------------------------------
    plotResilience(PCscoreTop, namesTop, datasetName, nIter, leaveOut, topN);
end


function PCscore = computePCscore(X, nPCs)
% Per-odour PC score for X (odours x fruits). Identical pipeline to the main
% analyses: within-fruit z-scoring -> PCA (fruits as observations) ->
% sum_{k=1..nPCs} |coeff(:,k)| * explained(k).

    % Within-fruit (row-wise) z-scoring: z-score down each fruit column.
    mu = mean(X, 1, 'omitnan');
    sd = std(X, 0, 1, 'omitnan');
    Xz = (X - mu) ./ sd;
    Xz(~isfinite(Xz)) = 0;            % guard constant/zero columns

    XX = Xz.';                        % fruits x odours (observations x variables)

    [coeff, ~, ~, ~, explained] = pca(XX);

    k = min(nPCs, size(coeff,2));     % fewer PCs may return for small subsets
    w = explained(1:k);
    PCscore = sum( abs(coeff(:,1:k)) .* w.', 2 );
end


function plotResilience(PCscore, names, datasetName, nIter, leaveOut, topN)
% Boxplot + jittered grey points, one point per iteration, top-N odours.

    nOdours = size(PCscore,1);

    Y = PCscore(:);
    G = repmat((1:nOdours)', size(PCscore,2), 1);

    figure('Color','w','Position',[100 100 max(900, 34*nOdours) 520]);

    greyBox   = [0.75 0.75 0.75];
    greyPoint = [0.45 0.45 0.45];

    boxplot(Y, G, 'Symbol','', 'Colors','k', 'Widths',0.6);
    hold on

    % Restyle: solid whiskers, no end-caps
    set(findobj(gca,'Tag','Upper Whisker'), 'LineStyle','-');
    set(findobj(gca,'Tag','Lower Whisker'), 'LineStyle','-');
    set(findobj(gca,'Tag','Upper Adjacent Value'), 'Visible','off');
    set(findobj(gca,'Tag','Lower Adjacent Value'), 'Visible','off');

    % Grey box fill with no edge stroke, then hide original box outlines
    hBoxes = findobj(gca,'Tag','Box');
    for i = 1:numel(hBoxes)
        patch(get(hBoxes(i),'XData'), get(hBoxes(i),'YData'), ...
              greyBox, 'FaceAlpha',0.5, 'EdgeColor','none');
    end
    set(hBoxes, 'Visible','off');
    uistack(findobj(gca,'Tag','Median'),'top');

    jitter = 0.18;
    for g = 1:nOdours
        yg = PCscore(g,:).';
        xg = g + (rand(numel(yg),1) - 0.5) * 2 * jitter;
        scatter(xg, yg, 14, ...
            'MarkerFaceColor', greyPoint, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.35);   % lower alpha: many points per odour
    end

    set(gca, 'XTick', 1:nOdours, 'XTickLabel', names, ...
             'XTickLabelRotation', 45, ...
             'TickLabelInterpreter','none', ...
             'TickDir','out', 'FontSize', 10, 'Box','off');
    xlim([0.5, nOdours+0.5]);
    ylabel('PC score');
    title(sprintf(['%s  —  leave-%d-fruits-out PC score  ' ...
                   '(top %d odorants; %d iterations)'], ...
                   datasetName, leaveOut, topN, nIter), 'Interpreter','none');
    hold off
end
