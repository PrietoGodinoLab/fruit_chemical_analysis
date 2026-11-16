%% ========================================================================
%  pcscore_dmel_no_citrus.m
%
%  CITRUS-EXCLUSION CONTROL 
%
% All citrus fruits
%  (grapefruit, mandarina, orange) are removed up front, and the PC score is
%  recomputed on the remaining fruits. 
%
%  The script does two things:
%    (1) Prints the rank and score of every limonene-containing compound on
%        the full panel and on the citrus-excluded panel (headline numbers
%        for the rebuttal letter).
%    (2) Runs the same leave-3-fruits-out resilience analysis as the main
%        Dmel script, but on the citrus-excluded panel, and plots it.
%
%  PC-score pipeline is identical to the main analyses (within-fruit z-scoring,
%  PCA with fruits as observations, sum of |loadings| * percent variance
%  explained over the first three PCs).
%  ========================================================================

clear; close all; clc;

%% ------------------------- USER SETTINGS --------------------------------
dmelFile     = 'dmel_data.csv';
citrusFruits = ["grapefruit","mandarina","orange"];  % columns to remove

nPCs        = 3;
leaveOut    = 3;
nIterWanted = 1000;
topN        = 25;
randSeed    = 1;

%% ------------------------- LOAD -----------------------------------------
T = readtable(dmelFile, 'VariableNamingRule','preserve');
odourNames = strtrim(string(T{:,1}));
fruitNames = string(T.Properties.VariableNames(2:end));
X          = T{:,2:end};
if ~isnumeric(X); X = str2double(string(X)); end

nFruits = size(X,2);

%% ------------------------- REMOVE CITRUS --------------------------------
isCitrus = ismember(lower(fruitNames), lower(citrusFruits));
if nnz(isCitrus) ~= numel(citrusFruits)
    warning(['Expected %d citrus columns but matched %d. ' ...
             'Check names. Matched: %s'], ...
            numel(citrusFruits), nnz(isCitrus), ...
            strjoin(fruitNames(isCitrus), ', '));
end
fprintf('Removing citrus fruits: %s\n', strjoin(fruitNames(isCitrus), ', '));

Xfull   = X;                      % 33 fruits
Xnoci   = X(:, ~isCitrus);        % citrus-excluded
fruitsNoci = fruitNames(~isCitrus);
fprintf('Full panel: %d fruits | Citrus-excluded: %d fruits\n', ...
        size(Xfull,2), size(Xnoci,2));

%% ------------------------- (1) LIMONENE RANK TABLE ----------------------
scFull = computePCscore(Xfull, nPCs);
scNoci = computePCscore(Xnoci, nPCs);

rankFull = tiedrank(-scFull);     % rank 1 = highest score
rankNoci = tiedrank(-scNoci);

limMask = contains(lower(odourNames), 'limonene');
limIdx  = find(limMask);

fprintf('\n--- Limonene-containing compounds: rank (score) ---\n');
fprintf('%-28s %14s %18s\n', 'compound', 'full panel', 'citrus-excluded');
for i = limIdx.'
    fprintf('%-28s   %3d (%6.2f)      %3d (%6.2f)\n', ...
        odourNames(i), rankFull(i), scFull(i), rankNoci(i), scNoci(i));
end

fprintf('\n--- Top 8 odorants, citrus-excluded panel ---\n');
[~, ord] = sort(scNoci, 'descend');
for r = 1:8
    i = ord(r);
    fprintf('  %2d. %7.2f  %s\n', r, scNoci(i), odourNames(i));
end

%% ------------------------- (2) LEAVE-3-OUT ON CITRUS-EXCLUDED -----------
runResilienceRandom(Xnoci, odourNames, 'D. melanogaster (citrus-excluded)', ...
                    nPCs, leaveOut, nIterWanted, topN, randSeed);


%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function runResilienceRandom(X, odourNames, datasetName, nPCs, leaveOut, ...
                             nIterWanted, topN, randSeed)
    nOdours = size(X,1);
    nFruits = size(X,2);

    rng(randSeed);
    maxCombos = nchoosek(nFruits, leaveOut);
    if nIterWanted >= maxCombos
        combos = nchoosek(1:nFruits, leaveOut);
        fprintf('\nEnumerating all %d leave-%d-out combinations.\n', maxCombos, leaveOut);
    else
        combos = zeros(nIterWanted, leaveOut);
        seen   = containers.Map('KeyType','char','ValueType','logical');
        count  = 0;
        while count < nIterWanted
            c   = sort(randperm(nFruits, leaveOut));
            key = sprintf('%d-', c);
            if ~isKey(seen, key)
                seen(key) = true; count = count + 1; combos(count,:) = c;
            end
        end
        fprintf('\nDrew %d unique random leave-%d-out combinations (of %d).\n', ...
                nIterWanted, leaveOut, maxCombos);
    end
    nIter = size(combos,1);

    PCscoreAll = nan(nOdours, nIter);
    for it = 1:nIter
        keep = true(1,nFruits); keep(combos(it,:)) = false;
        PCscoreAll(:,it) = computePCscore(X(:,keep), nPCs);
    end

    medScore   = median(PCscoreAll, 2, 'omitnan');
    [~, order] = sort(medScore, 'descend');
    keepIdx    = order(1:min(topN,nOdours));

    plotResilience(PCscoreAll(keepIdx,:), odourNames(keepIdx), ...
                   datasetName, nIter, leaveOut, topN);
end


function PCscore = computePCscore(X, nPCs)
    mu = mean(X, 1, 'omitnan');
    sd = std(X, 0, 1, 'omitnan');
    Xz = (X - mu) ./ sd;
    Xz(~isfinite(Xz)) = 0;
    XX = Xz.';
    [coeff, ~, ~, ~, explained] = pca(XX);
    k = min(nPCs, size(coeff,2));
    w = explained(1:k);
    PCscore = sum( abs(coeff(:,1:k)) .* w.', 2 );
end


function plotResilience(PCscore, names, datasetName, nIter, leaveOut, topN)
    nOdours = size(PCscore,1);
    Y = PCscore(:);
    G = repmat((1:nOdours)', size(PCscore,2), 1);

    figure('Color','w','Position',[100 100 max(900,34*nOdours) 520]);
    greyBox=[0.75 0.75 0.75]; greyPoint=[0.45 0.45 0.45];

    boxplot(Y, G, 'Symbol','', 'Colors','k', 'Widths',0.6); hold on

    % Restyle: solid whiskers, no end-caps
    set(findobj(gca,'Tag','Upper Whisker'), 'LineStyle','-');
    set(findobj(gca,'Tag','Lower Whisker'), 'LineStyle','-');
    set(findobj(gca,'Tag','Upper Adjacent Value'), 'Visible','off');
    set(findobj(gca,'Tag','Lower Adjacent Value'), 'Visible','off');

    % Grey box fill with no edge stroke, then hide original box outlines
    hBoxes = findobj(gca,'Tag','Box');
    for i = 1:numel(hBoxes)
        patch(get(hBoxes(i),'XData'), get(hBoxes(i),'YData'), ...
              greyBox,'FaceAlpha',0.5,'EdgeColor','none');
    end
    set(hBoxes, 'Visible','off');
    uistack(findobj(gca,'Tag','Median'),'top');

    jitter=0.18;
    for g=1:nOdours
        yg=PCscore(g,:).';
        xg=g+(rand(numel(yg),1)-0.5)*2*jitter;
        scatter(xg,yg,14,'MarkerFaceColor',greyPoint, ...
            'MarkerEdgeColor','none','MarkerFaceAlpha',0.35);
    end

    set(gca,'XTick',1:nOdours,'XTickLabel',names,'XTickLabelRotation',45, ...
        'TickLabelInterpreter','none','TickDir','out','FontSize',10,'Box','off');
    xlim([0.5,nOdours+0.5]); ylabel('PC score');
    title(sprintf(['%s  —  leave-%d-fruits-out PC score  ' ...
        '(top %d odorants; %d iterations)'], ...
        datasetName, leaveOut, topN, nIter),'Interpreter','none');
    hold off
end
