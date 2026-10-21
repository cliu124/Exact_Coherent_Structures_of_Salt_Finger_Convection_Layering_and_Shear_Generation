clear all;
close all;
clc;

set(0, 'Defaultfigurecolor',[1 1 1])
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 16)
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 16)

%% Load Reference CSV Data
ref_data1 = readmatrix('Holyer1984_Fig3.csv'); 
kz_ref1   = ref_data1(:, 1);
GR_ref1   = ref_data1(:, 2);

ref_data2 = readmatrix('Holyer1984_Fig6.csv'); 
kz_ref2   = ref_data2(:, 1);
GR_ref2   = ref_data2(:, 2);

%% Parameter Set 1 (Fig 3, Holyer1984)
Pr1 = 10.0; 
tau1 = 0.01; 
gamma1 = 0.5;
Rrho1 = gamma1 / tau1
Ra1 = gamma1 / (1 - gamma1)
aw0_1 = 4.0; 
max_kz1 = 0.5;

kz_list1 = linspace(0.001, max_kz1, 20)'; 
GR1 = getGrowthRate_kz(kz_list1, Ra1, Pr1, Rrho1, tau1, aw0_1);

% sample at kz = 0.3 for eigenvector visualization
kz_sample1 = 0.3;
[~,eigenval] = eig_Nguyen2026(aw0_1, Ra1, Pr1, Rrho1, tau1, kz_sample1);
GR1_sample = real(eigenval(1))

%% Parameter Set 2 (Fig 6, Holyer1984)
Pr2 = 1000.0; 
tau2 = 1.0 / 3.0; 
gamma2 = 0.91;
Rrho2 = gamma2 / tau2
Ra2 = gamma2 / (1 - gamma2)
aw0_2 = 100.0; 
max_kz2 = 0.25;

kz_list2 = linspace(0.001, max_kz2, 20)'; 
GR2 = getGrowthRate_kz(kz_list2, Ra2, Pr2, Rrho2, tau2, aw0_2);

%% Plotting Two-Panel Figure
% f = figure('Position', [100, 100, 1100, 450]);
% Set Units to 'inches' and Position to [x, y, width, height]
f = figure('Units', 'inches', 'Position', [1, 1, 6, 3]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% Panel (a)
ax1 = nexttile;
plot(kz_list1, real(GR1), 'k-', 'LineWidth', 2, 'DisplayName', 'Present'); hold on;
scatter(kz_ref1, GR_ref1, 50, 'r', 'o', 'filled', 'DisplayName', 'Holyer (1984)'); hold on;
scatter(kz_sample1, GR1_sample, 100, 'k', 'p', 'filled');
hold off;
xlabel('{\it{k_z}}');
ylabel('Growth rate');
xlim([0, max_kz1]);
% title('(a) \gamma = 0.5, Pr = 10, \tau = 0.01');
% legend('Location', 'northeast');

% Panel (b)
ax2 = nexttile;
plot(kz_list2, real(GR2), 'k-', 'LineWidth', 2, 'DisplayName', 'Present');
hold on;
scatter(kz_ref2, GR_ref2, 50, 'r', 'o', 'filled', 'DisplayName', 'Holyer (1984)');
hold off;
xlabel('{\it{k_z}}');
% ylabel('Growth rate');
xlim([0, max_kz2]);
% title('(b) \gamma = 0.91, Pr = 1000, \tau = 1/3');
legend('Location', 'east');

% Save entire figure
filename = 'validation_Holyer1984';
savefigure(f, [filename '.png']);

%% Save Function
function savefigure(fig_handle, filename)
    dpi = 300; 
    background_color = 'white'; 
    exportgraphics(fig_handle, filename, "Resolution", dpi, 'BackgroundColor', background_color)
end