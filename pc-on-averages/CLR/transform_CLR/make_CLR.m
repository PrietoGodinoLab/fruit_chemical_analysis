% --- Load original DMEL table (odours in col 1, samples in remaining cols)
T = readtable('dmel_data.csv');   % odour names in first column
odourNames = T{:,1};              % cell array / string array depending on file
sampleNames = T.Properties.VariableNames(2:end);

% Numeric matrix: rows = odours, cols = samples (as in your csv)
X = T{:,2:end};

% --- Step 1: set negative values to 0 (baseline artefacts)
X(X < 0) = 0;

% --- Step 2: transpose so rows = samples, cols = odours
X = X';   % now: samples x odours

% --- Step 3: zero replacement per sample
% Replace 0s with 0.5 * (smallest positive in that sample)
for i = 1:size(X,1)
    row = X(i,:);
    pos = row(row > 0);
    if ~isempty(pos)
        minPos = min(pos);
        row(row == 0) = 0.5 * minPos;
    else
        % If a sample somehow has no positives at all, leave as-is (rare)
        % (CLR would be undefined; you may want to remove this sample)
    end
    X(i,:) = row;
end

% --- Step 4: CLR transform
% CLR(x_j) = log(x_j / g(x)), where g(x) is geometric mean per sample
logX = log(X);
gm  = exp(mean(logX, 2));         % geometric mean per sample (samples x 1)
CLR = logX - log(gm);             % broadcasting (samples x odours)

% --- Step 5: write CLR matrix and save compound identity separately

nOdours = size(CLR,2);

% Create safe variable names: Odour_1, Odour_2, ...
safeNames = "Odour_" + string(1:nOdours);

% Build table
CLR_T = array2table(CLR, 'VariableNames', cellstr(safeNames));

% Add sample names
CLR_T = addvars(CLR_T, string(sampleNames(:)), ...
                'Before', 1, 'NewVariableNames', "Sample");

% Write CLR data
writetable(CLR_T, 'Dmel_CLR_transformed.csv');

% --- CRITICAL: save mapping key
NameMap = table(safeNames(:), string(odourNames(:)), ...
    'VariableNames', {'Variable','OriginalCompoundName'});

writetable(NameMap, 'Dmel_odour_name_key.csv');