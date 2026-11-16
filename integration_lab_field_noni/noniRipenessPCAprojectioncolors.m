%% this script is works specifically with the excel sheet provided - with 14 common odours of lab and field data- combined
% it performs a PCA of the lab data first, and using the mu (mean) and std
% (standard deviation) values from the lab data, the field data can be
% projected on the lab fruits data PC space. In this way we can find which
% are the variables (odours) that are imporant for the different ripeness
% stages of the lab data, and apply that to our field collected data. This
% can help determine the ripeness (according to chemical profile) of the
% field-collected data.


% import data from excel call it data
data = readtable('input.xlsx', 'Sheet', 'input_for_matlab');
tablet= data; % this contains field followed by lab samples
%% convert table to double, and extract lab  (auer) data and field data
auer= tablet{:,57:end}'; % this is to convert from table to matrix, and also exclude the column with odour names
field= tablet{:,2:56}';

%% calcualte relative abundance of each odor

 auer_per=(auer./sum(auer,2))*100;


%% z-norm auer (lab) data and perform PCA

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

auer_norm= zscor_xnan(auer_per')'; % zscore and puts in fruits vs odours table for pca

[coeff,score,latent,tsquared,explained,mu] = pca(auer_norm); % principal component analysis for n by p data matrix. Rows (n) corresponds to observations and columns (p) correspond to variables.
                 % coeff matrix is p by p. Each coumn of coeff contains
                 % coefficients for one pincipal component and the columns
                 % are in descending order of component variance. PCA
                 % default uses singular value decomposition algorithm. 

%% calculate the mean and std od z-normed lab (auer data). this will be used for PC projection

auer_mu= mean(auer_norm, 'omitnan');
auer_std= std(auer_norm,'omitnan');
                 
%% calculate variable coefficient weighted by variance explained
exp = explained'; %tranposes original explained vector, to match the PCs in the columns of coeff
coeffn= coeff.*exp; % matrix for coeffs weighted by variance explained

%% processing field data- normalise with auer data's mu and std.
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

field_per= (field./sum(field,2))*100; %relative percentage, abundance of odours within each fruit.
field_znorm= zscor_xnan(field_per')'; % z-normalise by fruit the field data

field_norm= (field_znorm- auer_mu )./ auer_std; % normalise with lab data's mu and std
field_project= field_norm*coeff; % projects normalised field data onto the PC space calculated for auer data
% by multiplying to the coefficient.

%% Plotting PCA scores in 3D of field data onto PC space of lab data (corrected data points)

% obtain average score for fruit repeats- to draw the line of best fit for
% ripening
avg_g= mean(score(1:3,:),'omitnan'); % green
avg_gy= mean(score(4:6,:), 'omitnan'); % green yellow
avg_y= mean(score(7:9,:));  % yellow
avg_yb= mean(score(10:12,:)); % yellow brown
avg_b= mean(score(13:15,:)); % brown

avg_auer= [avg_g; avg_gy; avg_y; avg_yb; avg_b]; % put mean scores into a vector

%3D plot starts here, dependency on rgb() function- to be downloaded 

figure
set(gcf, 'renderer', 'painters')

hold on
plot3(score(1:3,1), score(1:3,2),score(1:3,3),'.g','markersize',40)% G

plot3(score(4:6,1), score(4:6,2),score(4:6,3),'.','markersize',40, 'Color', [0.6, 0.80, 0.2]) % GY

plot3(score(7:9,1), score(7:9,2), score(7:9,3),'.y','markersize',40) % y

plot3(score(10:12,1), score(10:12,2), score(10:12,3),'.','markersize',40, 'Color', [0.61, 0.5, 0 ]) % YB

plot3(score(13:15,1), score(13:15,2),score(13:15,3),'.','markersize',40, 'Color', [0.58, 0.3, 0]) % B

plot3(field_project(1:5,1), field_project(1:5,2),field_project(1:5,3), '.', 'markersize', 40, 'color', [0.8235 0.7059 0.8392]) % s1

plot3(field_project([6:7 16:18],1), field_project([6:7 16:18],2),field_project([6:7 16:18],3), '.', 'markersize', 40, 'color', [0.4 0.8 0.67]) % s3

plot3(field_project([8 19:21],1), field_project([8 19:21],2), field_project([8 19:21],3), '.', 'markersize', 40, 'color', [0.545 0.529 0.271]) % s4

plot3(field_project([9:11 22:24],1), field_project([9:11 22:24],2),field_project([9:11 22:24],3), '.', 'markersize', 40, 'color', [1.0, 0.85, 0.73]) % s5

plot3(field_project(44:48,1), field_project(44:48,2),field_project(44:48,3), '.', 'markersize', 40, 'color', [0.5 0 1]) % s9 

plot3(field_project(49:52,1), field_project(49:52,2),field_project(49:52,3), '.', 'markersize', 40, 'color', [0.86, 0.08, 0.24]) % s10

plot3(field_project([25 53:55],1), field_project([25 53:55],2),field_project([25 53:55],3), '.', 'markersize', 40, 'color', [0.545, 0, 0.545]) % s11a

plot3(field_project([12:13 26:28],1), field_project([12:13 26:28],2),field_project([12:13 26:28],3), '.', 'markersize', 40, 'color', [0.2117, 0.9137, 0.7960] ) % s11b

plot3(field_project([14:15 29:31],1), field_project([14:15 29:31],2),field_project([14:15 29:31],3), '.', 'markersize', 40, 'color', [0 0 1]) % s12

plot3(field_project(32:35,1), field_project(32:35,2),field_project(32:35,3), '.', 'markersize', 40, 'color', [1 0.84 0]) % s14

plot3(field_project(36:39,1), field_project(36:39,2),field_project(36:39,3), '.', 'markersize', 40, 'color', [0.69 0.77 0.87]) % s15

plot3(field_project(40:43,1), field_project(40:43,2),field_project(40:43,3), '.', 'markersize', 40, 'color', [0.58, 0.44, 0.86]) % s19

fnplt(cscvn(avg_auer(:,1:3 )'),'k',1);% taking the average score of the each of the repeats, then using that data and 3 pcs, to plot the curve on top cannot make outputplot into dashed line


mult= 3; % multiplier for plotting coefficient vectors 

num= [4 6 7 8 9 10]; %index for the top PCA coeff odours
cl= {'SlateGray', 'darkseagreen', 'moccasin','RosyBrown', 'darkkhaki' ,'hotpink', 'lemonchiffon' }; % colours

%vector plots from origin for PCA coeffs of odours of interest
for i= 1:numel(num)
    plot3([0 coeff(num(i),1)].*mult, [0 coeff(num(i),2)].*mult, [0 coeff(num(i),3)], '-', 'Linewidth', 2, 'Color',rgb(cl(i)) ) % odours 
end 

hold off
grid on

xlabel('PC1- 67%') % titles- containing variance explained in each PC
ylabel('PC2- 21%')
zlabel('PC3- 11%')

legend('Green', 'Green-Yellow', 'Yellow', 'Yellow-Brown', 'Brown', 'Field 1', 'Field 3', 'Field 4', 'Field 5', 'Field 9' , 'Field 10', 'Field 11a', 'Field 11b', 'Field 12', 'Field 14', 'Field 15', ...
'Field 19', 'Methyl hexanoate', 'Ethyl hexanoate', 'Hexanoic acid', 'Methyl octanoate','Ethyl octanoate', 'Octanoic acid', 'Location', 'bestoutside') 



%% if you want to plot the loadings/coefficients of the PCA 

%% cummulative plot for variance explained in PCs

n= size(coeff,2); % no of PCs.
m = 1:n; % vector for number of PCs

figure
h = plot (cumsum(explained),'.b');
set(h,'MarkerSize',30);
set(h,'Clipping', 'off');
set(gcf,'color','w')
set(gca,'XTick', 1:n,'XTickLabel',m)
xlim([1,10])
ylim([0,100]);
set(gca,'TickDir','out')
set(gca,'FontSize',16)
set(gca,'box','off')
xlabel('PC')
ylabel('Cumulative variance explained')

%% PCA loadings * variance explained
figure

titl= "PCA of Lab noni" ;% change as data changes

k = stem(coeffn(:,1:size(coeffn,2))','.');
set(k,'MarkerSize', 30);
set(gcf,'color','w')
set(gca,'XTick', 1:3,'XTickLabel',m)
xlim([0.2,4.2])
set(gca,'TickDir','out')
set(gca,'FontSize',16)
set(gca,'box','off')
%legend ('show')
hold on
plot([0,100],[ 0 0], 'k-')
hold off
xlabel('PC')
ylabel('Loading*variance explained')
title(titl)

%% extract Odour names for the variables with highest loading*variance explained

odours= data{:,1}; % obtain the odours names from original table of data 
odcoeffn= table(odours, coeffn(:,1), coeffn(:,2), coeffn(:,3), coeffn(:,4)); % create table with odour names on the left side of the coeffs of each PC represented on the plot
% PC1
odmax_pc1= sortrows(odcoeffn,2); % pc1 sort min to max  

first_pc1= odmax_pc1{[1 end],1}; % gets the minimum and maximum values of pc1 (the odours accounting for highest variance)
second_pc1 = odmax_pc1{[2, end-1],1}; % gets SECOND min and max of pc1 (odours accounting second highest variance 
pc1= [first_pc1, second_pc1]  % displays first and second min max of PC1 set as:
                               % [ min PC1.   second min PC1
                                 % Max PC1    second max PC1] 

% PC2
odmax_pc2= sortrows(odcoeffn,3); % same as for PC1
first_pc2= odmax_pc2{[1 end],1};
second_pc2= odmax_pc2{[2 end-1],1};
pc2= [first_pc2, second_pc2]  % min and max in first and second pc. as well as first max/min and sec max/min 

% PC3 
odmax_pc3= sortrows(odcoeffn, 4);
first_pc3= odmax_pc3{[1 end],1};
second_pc3= odmax_pc3{[2, end-1],1};
pc3= [first_pc3, second_pc3]


%projection in PC1
figure
set(gcf, 'renderer', 'painters')

hold on
plot3(mean(score(1:3,1)), 0,0,'.g','markersize',40)% G

plot3(mean(score(4:6,1)), 0,0,'.','markersize',40, 'Color', [0.6, 0.80, 0.2]) % GY

plot3(mean(score(7:9,1)), 0, 0,'.y','markersize',40) % y

plot3(mean(score(10:12,1)), 0, 0,'.','markersize',40, 'Color', [0.61, 0.5, 0 ]) % YB

plot3(mean(score(13:15,1)), 0,0,'.','markersize',40, 'Color', [0.58, 0.3, 0]) % B

plot3(mean(field_project(1:5,1)), 0,0, '.', 'markersize', 40, 'color', [0.8235 0.7059 0.8392]) % s1

plot3(mean(field_project([6:7 16:18],1)), 0,0, '.', 'markersize', 40, 'color', [0.4 0.8 0.67]) % s3

plot3(mean(field_project([8 19:21],1)), 0, 0, '.', 'markersize', 40, 'color', [0.545 0.529 0.271]) % s4

plot3(mean(field_project([9:11 22:24],1)), 0,0, '.', 'markersize', 40, 'color', [1.0, 0.85, 0.73]) % s5

plot3(mean(field_project(44:48,1)), 0,0, '.', 'markersize', 40, 'color', [0.5 0 1]) % s9 

plot3(mean(field_project(49:52,1)), 0,0, '.', 'markersize', 40, 'color', [0.86, 0.08, 0.24]) % s10

plot3(mean(field_project([25 53:55],1)), 0,0, '.', 'markersize', 40, 'color', [0.545, 0, 0.545]) % s11a

plot3(mean(field_project([12:13 26:28],1)), 0,0, '.', 'markersize', 40, 'color', [0.2117, 0.9137, 0.7960] ) % s11b

plot3(mean(field_project([14:15 29:31],1)), 0,0, '.', 'markersize', 40, 'color', [0 0 1]) % s12

plot3(mean(field_project(32:35,1)), 0,0, '.', 'markersize', 40, 'color', [1 0.84 0]) % s14

plot3(mean(field_project(36:39,1)), 0,0, '.', 'markersize', 40, 'color', [0.69 0.77 0.87]) % s15

plot3(mean(field_project(40:43,1)), 0,0, '.', 'markersize', 40, 'color', [0.58, 0.44, 0.86]) % s19
