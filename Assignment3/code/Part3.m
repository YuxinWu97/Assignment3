%Part3
clear all 
clear

global C
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s
C.temp = 300;                       % Initial temperature 
C.m_e = 0.26*C.m_0;                 % Effective mass 


Wb = 5:5:50;
Lb = 40;

current = zeros(size(Wb,2)*size(Lb,2),1);
calcurrent = zeros(size(Wb,2)*size(Lb,2),1);
bnarea = zeros(size(Wb,2)*size(Lb,2),1);
                               

idx = 1;
for i= 1:size(Lb,2)
    for j=1:size(Wb,2)
        bnarea(idx,1) = Lb(i)*Wb(j);
        [current(idx,1), calcurrent(idx,1)] = den(Lb(i), Wb(j));
        idx = idx +1
    end
end


figure(1)
plot(Wb, calcurrent);
grid;
title('Current vs Bottleneck Width');
xlabel('Bottleneck Width (nm)');
ylabel('Current (A)');