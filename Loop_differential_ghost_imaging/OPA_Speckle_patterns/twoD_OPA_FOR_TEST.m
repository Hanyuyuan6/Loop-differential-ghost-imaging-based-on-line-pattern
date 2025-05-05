clc;    
clear;  

A = 1;         
Nx = 128;       
Ny = 128;     
lam = 1550e-9; 
Zx = 2000e-9;  
Zy = 2000e-9;  
a1 = 1000e-9;  
a2 = 1000e-9;   
k = 2 * pi / lam; 
neff = 1.56;    
nct = 1;        
thix = -4.05 * pi / 180 : 0.01 * pi / 180 : 4.05 * pi / 180; 
thiy = thix'; 

for i = 1:100
Eall = 0;  
for s = 1:Nx-1   
    E = 0;  
    E1 = exp(-1i * (s * (k .* Zx .* sin(thix) - (rand(1,1) - 0.5) .* 2 * pi))); 
    
    for h = 1:Ny-1 
        E2 = exp(-1i * (h * (k .* Zy .* sin(thiy) - (rand(1,1) - 0.5) .* 2 * pi)));
        E = E + E2;  
    end
    E = E1 .* E;  
    Eall = Eall + E;  
end

Eall = abs(Eall);  

T = A * Nx * Ny .* (sinc((pi .* a1) ./ lam .* sin(thix)) ...
    .* sinc((pi .* a2) ./ lam .* sin(thiy)) .* Eall) .^ 2;

thix0 = rad2deg(thix);
thiy0 = rad2deg(thiy);
figure(1); 
mesh(thix0, thiy0, T);  
xlabel('x direction (angle)'); 
ylabel('y direction (angle)');  
zlabel('Intensity (a.u.)');     
end
