clc;
clear;

% parameter definition
A = 1;
Nx = 128;
Ny = 128;
lam = 1550e-9;
Zx = 2000e-9;
Zy = 2000e-9;
a1 = 1000e-9;
a2 = 1000e-9;
k = 2*pi/lam;
neff = 1.56;
nct = 1;
thix = -4.05*pi/180 : 0.01*pi/180 : 4.05*pi/180;
thiy = thix';

for i = 1: 100
    lamx = 1000 * 10^(-9) * (rand(1,1) - 0.5) / 5 / 2;
    Eall = 0;
    for h = 0:Ny-1
        E = exp(1i * (k * h * Zy .* sin(thiy) - (rand(1,1) - 0.5) .* 2 * pi));
        Eall = Eall + E;
    end
    Eall = Eall .* sin(Nx/2 * (k * Zx .* sin(thix) - k * neff * lamx)) ./ ...
        sin(1/2 * (k * Zx .* sin(thix) - k * neff * lamx));
    Eall = abs(Eall);
    T = A * Nx * Ny .* ((sinc((pi * a1) ./ lam .* sin(thix))) .* ...
        (sinc((pi * a2) ./ lam .* sin(thiy))) .* Eall).^2;
    thix_0 = rad2deg(thix);
    thiy_0 = rad2deg(thiy);
    figure(1);
    mesh(thix_0, thiy_0, T);
    xlabel('x direction (angle)');
    ylabel('y direction (angle)');
    zlabel('Intensity (a.u.)');
end
