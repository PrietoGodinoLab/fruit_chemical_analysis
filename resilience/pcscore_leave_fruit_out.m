%% ========================================================================
%  pcscore_leave_fruit_out.m
%
%  INPUT-RESILIENCE ANALYSIS (leave-one-fruit-out PC score)
%
%  Purpose
%  -------
%  Tests whether the PC-score ranking of odorants is robust to which fruits
%  are included in the analysis. For a dataset of N fruits, the PCA + PC-score
%  pipeline is recomputed N times, each time leaving ONE fruit out. This yields
%  N PC-score values per odorant, whose spread shows how sensitive each
%  odorant's ranking is to the particular set of fruits sampled.
%
%  This is the leave-one-fruit-out analogue of the replicate-resampling
%  robustness analysis used elsewhere in the paper. The PC-score computation
%  itself is identical to that used for the main analyses (see
%  chem_pca_plus_PCscore.m): within-fruit (row-wise) z-scoring, PCA with fruits
%  as observations and odorants as variables, and PC score defined as the sum
%  of absolute PC loadings weighted by percent variance explained across the
%  first three principal components.
%
%  Input
%  -----
%  A CSV with odorants as ROWS and fruits as COLUMNS:
%     - column 1  : odorant names
%     - columns 2+: one column per fruit (already replicate-averaged)
%  Two datasets are analysed in turn: noni (D. sechellia) and pandan
%  (D. erecta). Each produces its own figure.
%
%  Output
%  ------
%  One figure per dataset: x-axis = odorants (sorted by median PC score,
%  descending), y-axis = PC score, one dot per leave-one-fruit-out iteration,
%  overlaid on a boxplot. All points shown in grey (no behavioural annotation).
%

%  ========================================================================

clear; close all; clc;

%% ------------------------- USER SETTINGS --------------------------------
% Edit these paths to point at your local copies of the two input files.
noniFile   = 'dsec_data_fruit_avgs.csv';
pandanFile = 'pandanfruit_averages.csv';

nPCs        = 3;      % number of leading PCs used in the PC score
leaveOut    = 1;      % number of fruits removed per iteration (1 = leave-one-out)
                      % (leaveOut > 1 handled by exhaustive combinations; see below)

%% ------------------------- RUN BOTH DATASETS ----------------------------
runResilience(noniFile,   'Noni (D. sechellia)',   nPCs, leaveOut);
runResilience(pandanFile, 'Pandan (D. erecta)',    nPCs, leaveOut);


%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function runResilience(csvFile, datasetName, nPCs, leaveOut)
% Load one dataset, run leave-fruit-out PC score, and plot.

    %% ---- Load data -----------------------------------------------------
    % Read with original column headers preserved; first column = odour names.
    T = readtable(csvFile, 'VariableNamingRule','preserve');

    odourNames = strtrim(string(T{:,1}));   % odorant names (rows)
    X          = T{:,2:end};                % odours x fruits (numeric matrix)

    % Coerce to double in case any column was read as non-numeric
    if ~isnumeric(X)
        X = str2double(string(X));
    end

    nOdours = size(X,1);
    nFruits = size(X,2);
    fruitNames = string(T.Properties.VariableNames(2:end));

    fprintf('\n%s: %d odorants x %d fruits\n', datasetName, nOdours, nFruits);

    %% ---- Enumerate leave-out iterations --------------------------------
    % Each row of "combos" lists the fruit indices REMOVED in that iteration.
    if leaveOut == 1
        combos = (1:nFruits)';                 % nFruits iterations, drop one each
    else
        combos = nchoosek(1:nFruits, leaveOut); % exhaustive combinations
    end
    nIter = size(combos,1);

    %% ---- Compute PC score for each iteration ---------------------------
    % PCscoreAll: odours x iterations
    PCscoreAll = nan(nOdours, nIter);
    for it = 1:nIter
        keep = true(1,nFruits);
        keep(combos(it,:)) = false;     % drop the selected fruit(s)
        Xsub = X(:, keep);              % odours x (nFruits - leaveOut)
        PCscoreAll(:,it) = computePCscore(Xsub, nPCs);
    end

    %% ---- Order odours by median PC score (descending) ------------------
    medScore = median(PCscoreAll, 2, 'omitnan');
    [~, order] = sort(medScore, 'descend');

    PCscoreOrdered = PCscoreAll(order, :);
    namesOrdered   = odourNames(order);

    %% ---- Plot ----------------------------------------------------------
    plotResilience(PCscoreOrdered, namesOrdered, datasetName, nIter);
end


function PCscore = computePCscore(X, nPCs)
% Compute per-odour PC score for a matrix X (odours x fruits).
% Reproduces the pipeline in chem_pca_plus_PCscore.m exactly:
%   1. within-fruit (row-wise) z-scoring of the odour profile
%   2. PCA with fruits as observations, odorants as variables
%   3. PCscore = sum over first nPCs of |loading| * percent-variance-explained
%
% Returns a column vector (odours x 1).

    % ---- Within-fruit z-scoring -------------------------------------
    % X is odours x fruits. We z-score DOWN each fruit column (i.e. across the
    % odours within that fruit), normalising each fruit's bouquet. NaNs ignored.
    % This matches: xx = zscor_xnan(x')'  with x = X'  in the original script.
    mu = mean(X, 1, 'omitnan');           % 1 x fruits
    sd = std(X, 0, 1, 'omitnan');         % 1 x fruits (N-1 normalisation)
    Xz = (X - mu) ./ sd;                  % odours x fruits, within-fruit z-scored

    % Guard against any all-constant fruit column (sd == 0 -> Inf/NaN)
    Xz(~isfinite(Xz)) = 0;

    % ---- Arrange for PCA: observations = fruits, variables = odours -----
    XX = Xz.';                            % fruits x odours

    % ---- PCA ------------------------------------------------------------
    % pca() centers columns (odours) internally, as in the original analysis.
    [coeff, ~, ~, ~, explained] = pca(XX);

    % Number of components actually returned can be < nPCs for small samples
    % (e.g. after leaving fruits out). Use whatever is available up to nPCs.
    k = min(nPCs, size(coeff,2));

    w = explained(1:k);                   % k x 1, percent variance explained
    PCscore = sum( abs(coeff(:,1:k)) .* w.', 2 );   % odours x 1
end


function plotResilience(PCscore, names, datasetName, nIter)
% Boxplot + jittered points of PC score per odour, one point per iteration.
% All grey (no behavioural annotation).
%
% PCscore : odours x iterations, already ordered for display
% names   : odours x 1 string, matching row order of PCscore

    nOdours = size(PCscore,1);

    % ---- Reshape to long form for grouped boxplot -----------------------
    % Y: all PC-score values stacked; G: group index (1..nOdours) per value
    Y = PCscore(:);                                  % (odours*iter) x 1
    G = repmat((1:nOdours)', nIter, 1);              % matching group indices

    figure('Color','w','Position',[100 100 max(900, 26*nOdours) 500]);

    greyBox   = [0.75 0.75 0.75];
    greyPoint = [0.45 0.45 0.45];

    % ---- Boxplot (no default outlier symbols; we draw all points) -------
    boxplot(Y, G, ...
        'Symbol','', ...
        'Colors','k', ...
        'Widths',0.6);
    hold on

    % ---- Restyle boxplot elements --------------------------------------
    % Whiskers: solid instead of dashed
    set(findobj(gca,'Tag','Upper Whisker'), 'LineStyle','-');
    set(findobj(gca,'Tag','Lower Whisker'), 'LineStyle','-');
    % Whisker end-caps: hide them
    set(findobj(gca,'Tag','Upper Adjacent Value'), 'Visible','off');
    set(findobj(gca,'Tag','Lower Adjacent Value'), 'Visible','off');

    % Colour the box interiors grey (no edge stroke), then hide the
    % original black box outlines so only the grey fill remains.
    hBoxes = findobj(gca,'Tag','Box');
    for i = 1:numel(hBoxes)
        patch(get(hBoxes(i),'XData'), get(hBoxes(i),'YData'), ...
              greyBox, 'FaceAlpha',0.5, 'EdgeColor','none');
    end
    set(hBoxes, 'Visible','off');

    % Redraw median line on top of the patches
    uistack(findobj(gca,'Tag','Median'),'top');

    % ---- Jittered individual points (one per iteration) -----------------
    jitter = 0.18;
    for g = 1:nOdours
        yg = PCscore(g,:).';
        xg = g + (rand(numel(yg),1) - 0.5) * 2 * jitter;
        scatter(xg, yg, 22, ...
            'MarkerFaceColor', greyPoint, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.6);
    end

    % ---- Cosmetics ------------------------------------------------------
    set(gca, 'XTick', 1:nOdours, 'XTickLabel', names, ...
             'XTickLabelRotation', 45, ...
             'TickLabelInterpreter','none', ...
             'TickDir','out', 'FontSize', 10, 'Box','off');
    xlim([0.5, nOdours+0.5]);
    ylabel('PC score');
    title(sprintf('%s  —  leave-one-fruit-out PC score (%d iterations)', ...
                  datasetName, nIter), 'Interpreter','none');
    hold off
end
