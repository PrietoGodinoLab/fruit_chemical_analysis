%% PC_score_robustness script 
% This script takes a random combination of sample repeats from each fruit
% to calculate the PCA (and subsequently outputs PC coefficients/loadings 
% of the first 3PCs) of that set of samples.
% User can select how many unique combinations to run. 

% INPUT: gc-data that is already thresholded (>10000), with odours in rows,
% and samples in columns (containing all sample repeats and clearly
% labelled sample names)

% OUTPUT: PC coefficients/loadings table, with odours in rows, for each unique combination
% of fruit samples selected in columns. Adds the magnitudes of the coeffs
% of the first 3 PCs.
% Also outputs a box plot plotting the above data.


%% Load data and set the output file names (variables to change before running the script)
data = readtable(['input_Dsec_robustness.csv']);

%set the file names for outputs
filename_pcscores= "output_PCscore_noni.csv";
filename_combos = "output_PCscore_robustness_combos_noni_metadata.csv";

plot_name = 'PCscore_robustness_boxplot_pandan.pdf';
fruit_name= 'Pandan'; % used in the plot title

% number combos wanted:
numSamples = 1000;



%% Group columns by fruit
colnames = data.Properties.VariableNames(2:end); % skip Odour
fruits = cellfun(@(x) regexp(x,'^([^_]+)','tokens','once'), colnames);
%fruits = cellfun(@(x) x{1}, fruits, 'UniformOutput', false);

uniqueFruits = unique(fruits);
fruitGroups = containers.Map;
for i = 1:numel(uniqueFruits)
    fruit = uniqueFruits{i};
    reps = colnames(strcmp(fruits, fruit));
    fruitGroups(fruit) = reps;
end
%% Generate random combinations (with seed)
rng(100); % reproducibility

randomCombos = cell(numSamples,1);
for k = 1:numSamples
    combo = cell(size(uniqueFruits));
    for i = 1:numel(uniqueFruits)
        reps = fruitGroups(uniqueFruits{i});
        r = randi(numel(reps)); % random replicate
        combo{i} = reps{r};
    end
    randomCombos{k} = combo;
end

%% Remove duplicate combinations before PCA

% Convert each combo (cell array of sample names) into a single string
comboStrings = cellfun(@(c) strjoin(c, ','), randomCombos, 'UniformOutput', false);

% Find unique combinations
[uniqueComboStrings, ia] = unique(comboStrings, 'stable');
uniqueCombos = randomCombos(ia);

% Report uniqueness stats
numTotal  = numel(comboStrings);
numUnique = numel(uniqueComboStrings);
fprintf('Total combos generated: %d\n', numTotal);
fprintf('Unique combos retained: %d\n', numUnique);
fprintf('Removed %d duplicate combos.\n', numTotal - numUnique);

% Update number of samples
numSamples = numUnique;
randomCombos = uniqueCombos;


%% Run PCA for each random combo and collect results
pcscore_wide = table(data.Odour, 'VariableNames', {'Odour'});
combo_metadata = cell2table(cell(numSamples, numel(uniqueFruits)), ...
    'VariableNames', uniqueFruits);

for i = 1:numSamples
    combo = randomCombos{i};
    selectedCols = [{'Odour'}, combo];
    df_combo = data(:, selectedCols);
    
    % Store metadata
    combo_metadata{i,:} = combo;
    
    % Run PCA + PC score
    pc_df = run_pca_and_pcscore(df_combo);
    pcscore_wide.(['Combo_' num2str(i)]) = pc_df.PC_score;
end


%% Save results
writetable(pcscore_wide, filename_pcscores);
writetable(combo_metadata, filename_combos);


%% Summary boxplot (sorted by max PC score, descending)

% Numeric matrix: rows = odours, cols = combos
X_raw = table2array(pcscore_wide(:,2:end)); % size: nOdours x nCombos

% Compute max score per odour (across combos)
%max_scores = max(X_raw, [], 2);

%compute median score per odour (across combos)
max_scores = median(X_raw, 2);

% Sort indices by descending max score
[~, sortIdx] = sort(max_scores, 'descend');

% Reorder data and labels
X_sorted = X_raw(sortIdx, :).';         % transpose so columns = odours
labels_sorted = pcscore_wide.Odour(sortIdx);

% Ensure labels are cell array of char
if isstring(labels_sorted)
    labels_sorted = cellstr(labels_sorted);
elseif ~iscell(labels_sorted)
    labels_sorted = cellstr(string(labels_sorted));
end

% Create figure
fig=figure('Position',[100 100 1600 900]);

% Boxplot with sorted odours
boxplot(X_sorted, 'Labels', labels_sorted, 'PlotStyle', 'compact');

% Formatting
xtickangle(90);
set(gca, 'FontSize', 8);
set(gca, 'TickLabelInterpreter', 'none');
xlabel('Odour (sorted by max PC score)');
ylabel('PC Score (|coeffn| sum of first 3 PCs)');
title(['Variation of PC Scores across Random Replicate Selections (',num2str(numSamples),') -', fruit_name])

% save the plot
exportgraphics(fig, plot_name, 'ContentType', 'vector');

%% Function to run PCA and compute PC score
function pc_df = run_pca_and_pcscore(data)
    x= data{:,2:end}';
    zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

    xx= zscor_xnan(x')'; % does the following steps in one line basically does the zscore and puts in fruits vs odours table for pca
    % a = x';
    % aa = zscor_xnan (a);
    %    xx = aa';
    %xx2= zscor_xnan(x2')'; 
    [coeff,score,latent,tsquared,explained,mu] = pca(xx); % principal component analysis for n by p data matrix. Rows (n) corresponds to observations and columns (p) correspond to variables.
                 % coeff matrix is p by p. Each coumn of coeff contains
                 % coefficients for one pincipal component and the columns
                 % are in descending order of component variance. PCA
                 % default uses singular value decomposition algorithm. 
                 
    
    % coeffn = coeff * explained
    coeffn = zeros(size(coeff));
    for i = 1:size(coeff,2)
        coeffn(:,i) = coeff(:,i) * explained(i);
    end
    
    % PC score = sum of abs(coeffn across first 3 PCs)
    pc_score = sum(abs(coeffn(:,1:3)),2);
    
    pc_df = table(data.Odour, pc_score, 'VariableNames', {'Odour','PC_score'});
end