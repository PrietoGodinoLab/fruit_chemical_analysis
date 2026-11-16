%% Chemical profiles: within-fruit z-scoring, PCA and odour PC scores
% Input CSV format: column 1 = odour names; remaining columns = fruit samples.
% Each row corresponds to an odour, each sample column to a fruit profile.
%
% Requires MATLAB Statistics and Machine Learning Toolbox (pca, boxplot).
% Behavioural annotations and the behavioural boxplot are used ONLY when
% inputFile is named exactly 'dmel_data.csv'. All other datasets run without
% the annotation workbook.
%
% Change this one line to select the dataset:
%inputFile = 'pandanfruit_averages.csv';
% inputFile = 'dsec_data_fruit_avgs.csv';
 inputFile = 'dmel_data.csv';

behaviourFile = 'behaviour-anotations.xlsx'; % used only for dmel_data.csv
showSampleLabels = true;                    % labels on PC scatterplots

%% 1. Read input data
if ~isfile(inputFile)
    error('Input file not found: %s', inputFile);
end

tablet = readtable(inputFile);
odourNames = strtrim(string(tablet{:,1}));
sampleNames = string(tablet.Properties.VariableNames(2:end));

% Transpose: rows = fruit samples, columns = odours (as required for PCA).
x = tablet{:,2:end}.';
if ~isnumeric(x)
    error('The abundance columns in %s must contain numeric data.', inputFile);
end
if size(x,1) < 3 || size(x,2) < 3
    error('At least three fruit samples and three odours are required.');
end

%% 2. Within-fruit (row-wise) z-scoring
% Each fruit is centred and scaled across its odours, not across fruits.
% This preserves the calculation in the original analysis script.
fruitMean = mean(x, 2, 'omitnan');
fruitStd  = std(x, 0, 2, 'omitnan');
if any(~isfinite(fruitStd) | fruitStd == 0)
    error('At least one fruit has zero/undefined within-fruit standard deviation.');
end
xx = bsxfun(@rdivide, bsxfun(@minus, x, fruitMean), fruitStd);

%% 3. PCA: fruits = observations, odours = variables
[coeff, score, latent, tsquared, explained, mu] = pca(xx);
nPC = size(coeff,2);
nScorePC = min(3, nPC);
if nScorePC < 3
    warning('Only %d PCs are available; PC scores use all available PCs.', nScorePC);
end

%% 4. PC score: sum of absolute loadings weighted by explained variance
% j = odour (row of coeff), i = principal component (column of coeff).
% 'explained' is already expressed in percentage points, as in the paper.
coeffn = bsxfun(@times, coeff, explained(:).');
PCscore = sum(abs(coeffn(:,1:nScorePC)), 2);

PCscoreTable = table(odourNames, PCscore, ...
    'VariableNames', {'Odour','PC_score_rowZ'});

fprintf('\nDataset: %s\n%d fruit samples; %d odours.\n', ...
    inputFile, size(x,1), numel(odourNames));
fprintf('Variance explained by first %d PC(s): %.2f%%\n', ...
    nScorePC, sum(explained(1:nScorePC)));

%% 5. Behavioural annotation plot: D. melanogaster ONLY
% This check uses the input filename (without its folder path).
[~, inputName, inputExtension] = fileparts(inputFile);
isDmel = strcmpi([inputName inputExtension], 'dmel_data.csv');

if isDmel
    if ~isfile(behaviourFile)
        error('Behavioural annotation file not found: %s', behaviourFile);
    end
    ann = readtable(behaviourFile);
    if ~all(ismember({'Odour','behaviour'}, ann.Properties.VariableNames))
        error('Annotation file must contain columns named Odour and behaviour.');
    end
    ann.Odour = strtrim(string(ann.Odour));

    % Match compounds by name; use only 1 (attractive), 0 (neutral),
    % and 0.5 (conflicting). Exclude -1 (aversive) and unannotated odours.
    M = innerjoin(ann, PCscoreTable, 'Keys', 'Odour');
    M = M(ismember(M.behaviour, [1 0 0.5]), :);
    if isempty(M)
        warning('No recognised behavioural annotations matched the input.');
    else
        behLabel = strings(height(M),1);
        behLabel(M.behaviour == 1)   = "Attractive";
        behLabel(M.behaviour == 0)   = "Neutral";
        behLabel(M.behaviour == 0.5) = "Conflicting";
        behLabel = categorical(behLabel, ...
            ["Attractive","Neutral","Conflicting"]);

        figure('Color','w', 'Name','Behavioural annotations');
        boxplot(M.PC_score_rowZ, behLabel, 'Symbol','');
        hold on;
        % Deterministic jitter for reproducible graphics.
        rng(1, 'twister');
        xpos = grp2idx(behLabel);
        jitterAmount = 0.18;
        xjitter = xpos + (rand(size(xpos)) - 0.5) * 2 * jitterAmount;
        scatter(xjitter, M.PC_score_rowZ, 35, 'o', 'filled', ...
            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','k');
        ylabel('Within-fruit z-score PC score (PC1-PC3)');
        title('Behaviour vs PC score');
        ylim([0 19]); 
        box on;
        set(gca,'FontSize',12);
        hold off;
        fprintf('Matched annotated odours: %d\n', height(M));
    end
else
    fprintf('Behavioural annotation plot skipped (not dmel_data.csv).\n');
end

%% 6. Cumulative variance explained
figure('Color','w', 'Name','Cumulative explained variance');
plot(1:nPC, cumsum(explained), '.b', 'MarkerSize',20);
xlim([0.5, min(10,nPC)+0.5]);
xticks(1:min(10,nPC));
ylim([0 100]);
xlabel('PC');
ylabel('Cumulative variance explained (%)');
box off;
set(gca,'FontSize',14, 'TickDir','out');

%% 7. Explained-variance-weighted loadings (original plot style)

nPlotPC = min(3, nPC);
figure('Color','w', 'Name','Variance-weighted loadings');
k = stem(coeffn(:,1:nPlotPC).', '.');
set(k, 'MarkerSize', 30);
set(gca, 'XTick', 1:nPlotPC, 'XTickLabel', 1:nPlotPC);
xlim([0.2, 3.2]);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 16);
set(gca, 'Box', 'off');
hold on;
plot([0, 100], [0, 0], 'k-');
hold off;
xlabel('PC');
ylabel('Loading*variance explained');

%% 8. Most positive and most negative loadings for PC1-PC3
% The sign of a PCA axis is arbitrary; both extremes are shown.
for pc = 1:nScorePC
    [~, order] = sort(coeff(:,pc), 'descend');
    fprintf('\nPC%d: most positive loadings\n', pc);
    for rank = 1:min(2,numel(order))
        idx = order(rank);
        fprintf('  %s: %.5g\n', char(odourNames(idx)), coeff(idx,pc));
    end
    fprintf('PC%d: most negative loadings\n', pc);
    for rank = 1:min(2,numel(order))
        idx = order(end-rank+1);
        fprintf('  %s: %.5g\n', char(odourNames(idx)), coeff(idx,pc));
    end
end

%% 9. PC scatterplots: use actual number and names of samples

figure('Color','w', 'Name','PC1 vs PC2');
scatter(score(:,1), score(:,2), 45, 'filled');
xlabel('PC1');
ylabel('PC2');
box on;
if showSampleLabels
    text(score(:,1), score(:,2), cellstr(" " + sampleNames(:)), ...
        'FontSize',8, 'Interpreter','none');
end

if nPC >= 3
    figure('Color','w', 'Name','PC1 vs PC2 vs PC3');
    scatter3(score(:,1), score(:,2), score(:,3), 45, 'filled');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    grid on;
    if showSampleLabels
        text(score(:,1), score(:,2), score(:,3), ...
            cellstr(" " + sampleNames(:)), ...
            'FontSize',8, 'Interpreter','none');
    end
end
