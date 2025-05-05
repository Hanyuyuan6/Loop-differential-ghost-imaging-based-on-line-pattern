clear;            
clc;              

lam = 1550e-9;    
a = 1000e-9;     
d = 2000e-9;     
N = 128;          
k = 2*pi/lam;     % k = 2π/λ
p = 0;            
phi = -4.05*pi/180 : 0.01*pi/180 : 4.05*pi/180;  
Rphi = [];        

for i = 1:100
    for s = 0:N-1
        % rphi = (rand(1,1)-0.5);
        rphi = randn(1, 1, 1);
        p = p + exp(1i*(k*s*d.*sin(phi) - rphi.*2*pi));
        i = s + 1;
        Rphi(i) = -rphi/2;
    end
    p = abs(p);
    P = N^2 .* (sinc((pi.*a)./lam.*sin(phi)).^2) .* p.^2;

    angle = rad2deg(phi);
    figure(1);
    plot(angle, P);
    xlabel('Position');
    ylabel('Intensity');
end
aaaa = P';
p = p';