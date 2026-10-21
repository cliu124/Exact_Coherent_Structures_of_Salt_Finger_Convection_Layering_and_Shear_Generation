function [eigenvector,eigenvalue] = eig_Nguyen2026(aw0,Ra,Pr,Rrho,tau,kz,N)
    arguments
        aw0; Ra; Pr; Rrho; tau; kz;
        N = 128; % Sets 128 as the default if N is omitted
    end
    % This function computes eigenvalues of linear equations of salt-finger convection 
    % in the vertically periodic domain.
    %
    % Define input paramesters:
    % aw0 - base flow amplitude
    % Ra - Rayleigh number
    % Pr - Prandtl number
    % Rrho - density ratio
    % tau - diffusivity ratio
    % kx, ky - horizontal wavenumbers
    % N - Fourier modes
    % Define output paramesters:
    % eigenvector and eigenvalues
    % 
    % N = 128; % number of Fourier modes
    if kz == 0
        error('Vertical wavenumber kz cannot be zero for eig_Nguyen2026 function.');
    elseif kz < 1e-3
        N = 256; % increase Fourier modes for small kz to ensure accuracy
    elseif kz < 1e-4
        N = 512;
    end

    n_sf = 1; % number of salt fingers in the horizontal domain
    % critical wavenumber for steady elevator mode
    kxc = (Ra*(1.0/(tau*Rrho)-1.0))^(1.0/4.0);
    lxc = n_sf*2*pi/kxc;
    
    % Identity and zero matrices
    I=eye(N,N);
    O=zeros(N,N);

    % Differentiation Matrices 
    [~, D1] = fourdif(N, 1);
    [x, D2] = fourdif(N, 2);

    % Coordinate transformation
    x = x/(kxc/n_sf); % rescale x in [0, 2*pi) to be in [0, lxc)
    D1 = (kxc/n_sf)*D1;
    D2 = (kxc/n_sf)^2*D2;

    % Differential operators
    DX = D1;
    DZ = 1i*kz*I;
    Lap = -kz^2*I + D2; % Laplacian

    % Base flow operators
    at0 = aw0/(kxc^2);
    as0 = aw0/(tau*kxc^2);
    wb = -aw0*sin(kxc*x);
    Dwb = -aw0*kxc*cos(kxc*x);
    tb = at0*sin(kxc*x);
    Dtb = at0*kxc*cos(kxc*x);
    sb = as0*sin(kxc*x);
    Dsb = as0*kxc*cos(kxc*x);

    % Construct A matrix (Right hand side: dX/dt = AX)
    A = [diag(-wb)*DZ+Pr*Lap, O                  , -DX, O               , O                   ;%u
         diag(-Dwb)         , diag(-wb)*DZ+Pr*Lap, -DZ, Pr*Ra*I         , -Pr*Ra/Rrho*I       ;%w
         DX                 , DZ                 , O  , O               , O                   ;%p
         diag(-Dtb)         , -I                 , O  , -diag(wb)*DZ+Lap, O                   ;%T
         diag(-Dsb)         , -I                 , O  , O               ,-diag(wb)*DZ+tau*Lap];%S
    B = blkdiag(I, I, O, I, I); % u w p T S

    [eigenvector,eigenvalue] = eig(A,B);
    eigenvalue = diag(eigenvalue);
    clean_mask = ~isinf(eigenvalue) & ~isnan(eigenvalue) & (abs(eigenvalue) < 1e2);
    eigenvalue = eigenvalue(clean_mask);
    eigenvector = eigenvector(:, clean_mask);
    [~, sort_indices] = sort(real(eigenvalue), 'descend');
    eigenvalue = eigenvalue(sort_indices);
    eigenvector = eigenvector(:, sort_indices);
end