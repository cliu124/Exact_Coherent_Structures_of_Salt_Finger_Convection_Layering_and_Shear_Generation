clear all;
close all;
clc;

set(0, 'Defaultfigurecolor',[1 1 1])
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 16)
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 16)

% Generate a 256-color Blue-White-Red colormap
n = 256;
half = n / 2;
% Blue to White, then White to Red
bwr_cmap = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1); ...
            ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];


% Define 11 ColorBrewer control points for PuOr_r (Purple -> White -> Orange)
nodes_PuOr_r = [
    0.176  0.000  0.294; % Dark Purple
    0.329  0.153  0.533;
    0.502  0.451  0.675;
    0.698  0.671  0.824;
    0.847  0.855  0.922;
    0.969  0.969  0.969; % Center White
    0.996  0.878  0.714;
    0.992  0.722  0.388;
    0.878  0.510  0.078;
    0.702  0.345  0.024;
    0.498  0.231  0.031  % Dark Orange
];
% Interpolate to 256 levels
N = 256;
PuOr_r = interp1(linspace(0, 1, 11), nodes_PuOr_r, linspace(0, 1, N));

%% Parameter Set 1 (Fig 3, Holyer1984)
Pr = 10.0; 
tau = 0.01; 
gamma = 0.5;
Rrho = gamma / tau;
Ra = gamma / (1 - gamma);
aw0 = 4.0; 
kz = 0.3;

N = 128; % number of Fourier modes

[eigenvec, eigenval] = eig_Nguyen2026(aw0, Ra, Pr, Rrho, tau, kz, N);

fprintf(['Eigenvalue for kz = ' num2str(kz, '%.2f') ': ' num2str(real(eigenval(1)), '%.6f') ' + ' num2str(imag(eigenval(1)), '%.6f') 'i\n']);

eigenvector = eigenvec(:, 1); % Select the eigenvector corresponding to the first eigenvalue
fprintf('Size of eigenvector: %d x %d\n', size(eigenvector));


n_sf = 1; % number of salt fingers in the horizontal domain
kx = (Ra*(1.0/(tau*Rrho)-1.0))^(1.0/4.0); % critical wavenumber for steady elevator mode
lxc = n_sf*2*pi/kx;

x=linspace(0.,lxc,N); % wavelength = 2*pi/kx
z=linspace(0.,2*pi/kz,N);

% eigenvector is organized as [u; w; p; T; S], each of size N
u = zeros(N,N);
w = zeros(N,N);
T = zeros(N,N);
S = zeros(N,N);
uvec = eigenvector(0*N+1:0*N+N);
wvec = eigenvector(1*N+1:1*N+N);
Tvec = eigenvector(3*N+1:3*N+N);
Svec = eigenvector(4*N+1:4*N+N);

as0 = aw0/(tau*kx^2);
sb = as0*sin(kx*x);

base_S = zeros(N,N);
total_S = zeros(N,N);

for x_index=1:N
    for z_index=1:N
        x_val = x(x_index);
        z_val = z(z_index);
        u(x_index,z_index) = real(exp(1i*kz*z_val)*uvec(x_index));
        w(x_index,z_index) = real(exp(1i*kz*z_val)*wvec(x_index));
        T(x_index,z_index) = real(exp(1i*kz*z_val)*Tvec(x_index));
        S(x_index,z_index) = real(exp(1i*kz*z_val)*Svec(x_index));

        base_S(x_index,z_index) = sb(x_index);
        total_S(x_index,z_index) = sb(x_index) + 100*real(exp(1i*kz*z_val)*Svec(x_index));
    end
end
maxu = max(u,[],"all");
maxw = max(w,[],"all");
maxT = max(T,[],"all");
maxS = max(S,[],"all");

%% plot 
f1 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolor(x,z,u'); shading interp;
colormap(PuOr_r);
colorbar;
title('$\hat{u}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_u.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
pcolor(x,z,w'); shading interp;
colormap(PuOr_r);
colorbar;
title('$\hat{w}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_w.png']);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
pcolor(x,z,T'); shading interp;
colormap(PuOr_r);
colorbar;
title('$\hat{T}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_T.png']);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
pcolor(x,z,S'); shading interp;
colormap(PuOr_r);
colorbar;
title('$e_1$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_S.png']);

% Save total salinity field
total_S = base_S + 100*S;
f5 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
pcolor(x,z,total_S'); shading interp;
colormap(PuOr_r);
colorbar;
title('$s+100e_1$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_total_S.png']);

f6 = figure('Units', 'inches', 'Position', [1, 1, 4, 3]);
pcolor(x,z,base_S'); shading interp;
colormap(PuOr_r);
colorbar;
title('$s$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk_x}'})
yticks([0 2*pi/kz])
yticklabels({'0','2\pi/{\itk_z}'})
daspect([1 1 1])
savefigure(gca,['validation_Holyer1984_eigenvector_base_S.png']);

%% Save Function
function savefigure(fig_handle, filename)
    dpi = 300; 
    background_color = 'white'; 
    exportgraphics(fig_handle, filename, "Resolution", dpi, 'BackgroundColor', background_color)
end