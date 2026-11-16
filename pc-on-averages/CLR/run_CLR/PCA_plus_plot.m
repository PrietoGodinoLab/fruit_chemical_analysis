%% Reproduce CLR PCA + behaviour boxplot
clear; clc;

%% ---------- Inputs ----------
clrFile  = 'Dmel_CLR_transformed.csv';
keyFile  = 'Dmel_odour_name_key.csv';
annFile  = 'behaviour-anotations.xlsx';

%% ---------- Load CLR matrix ----------
CLR_table = readtable(clrFile);

% rows = samples, cols = odours (Odour_1..Odour_N)
CLR = CLR_table{:, 2:end};

%% ---------- Load odour key (Odour_# -> real compound name) ----------
key = readtable(keyFile);
odourNames = string(key.OriginalCompoundName);  % length should equal #odours in CLR

%% ---------- PCA on CLR (MATLAB default SVD PCA) ----------
[coeff, score, latent, tsquared, explained, mu] = pca(CLR);

fprintf('CLR PCA variance explained:\n');
fprintf('  PC1: %.4f%%\n', explained(1));
fprintf('  PC2: %.4f%%\n', explained(2));
fprintf('  PC3: %.4f%%\n', explained(3));

%% ---------- Abs variance-weighted PC score (first 3 PCs) ----------
weights = explained(1:3);  % percent variance
Abs_PC_score_CLR = sum(abs(coeff(:,1:3)) .* weights.', 2);  % odours x 1

% Table of CLR scores keyed by real odour names
CLRscores = table(odourNames, Abs_PC_score_CLR, ...
    'VariableNames', {'Odour','Abs_PC_score_CLR'});

%% ---------- Load behavioural annotations ----------
ann = readtable(annFile);  % expects columns: Odour, PC_score, behaviour

% Basic cleanup (helps if there are trailing spaces etc.)
ann.Odour = string(ann.Odour);
CLRscores.Odour = string(CLRscores.Odour);

ann.Odour = strtrim(ann.Odour);
CLRscores.Odour = strtrim(CLRscores.Odour);

%% ---------- Merge annotations with CLR scores ----------
M = innerjoin(ann, CLRscores, 'Keys', 'Odour');

% Remove missing behaviour labels (if any)
M = M(~isnan(M.behaviour), :);

% Ignore aversive (-1) as you requested
M = M(M.behaviour ~= -1, :);

%% ---------- Create category labels ----------
% behaviour codes: 1 attractive, 0 neutral, 0.5 conflicting
behLabel = strings(height(M),1);
behLabel(M.behaviour == 1)   = "Attractive";
behLabel(M.behaviour == 0)   = "Neutral";
behLabel(M.behaviour == 0.5) = "Conflicting";

% Keep only rows we successfully labeled
keep = behLabel ~= "";
M = M(keep,:);
behLabel = categorical(behLabel(keep), ["Attractive","Neutral","Conflicting"]);

%% ---------- CLR Boxplot + individual datapoints ----------
figure('Color','w');

% Draw boxplot
boxplot(M.Abs_PC_score_CLR, behLabel, 'Symbol',''); % remove default outlier stars
hold on

% Convert category labels to numeric positions (1,2,3)
xpos = grp2idx(behLabel);

% Horizontal jitter so points don't overlap
jitterAmount = 0.18;
xjitter = xpos + (rand(size(xpos)) - 0.5) * 2 * jitterAmount;

% Plot each odour as a circle
scatter(xjitter, M.Abs_PC_score_CLR, ...
        35, ...                % marker size
        'o', ...
        'filled', ...
        'MarkerFaceAlpha',0.6,...
        'MarkerEdgeColor','k');

% Labels
ylabel('CLR absolute PC score (PC1–PC3)');
title('Behaviour vs CLR PC score');

% IMPORTANT: same scale as row-zscore plot
ylim([0 19]);

box on
set(gca,'FontSize',12)

hold off

%% ---------- Optional: export the merged table ----------
% This is the “similar table” you asked for: Odour | CLR PC score | behaviour
Out = M(:, {'Odour','Abs_PC_score_CLR','behaviour'});
Out = sortrows(Out, 'Abs_PC_score_CLR', 'descend');
writetable(Out, 'Annotated_CLR_PC_scores_from_MATLAB.csv');

disp('Top 10 annotated odours by CLR PC score:');
disp(Out(1:min(10,height(Out)),:));