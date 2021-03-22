%Part1
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
Econ = 1e15;
Econ = Econ/0.001;


xmax = 200e-9;
ymax = 100e-9;

tStep = .01e-12;

nStep = 150;
nAtom = 3000;
nPlot = 25;

dens = zeros(nStep,1);
drift = zeros(nStep,1);
time = linspace(0,tStep*nStep, nStep);

Vx = 0.1;
Vy = 0.1; 

Ex = Vx/xmax;
Ey = Vy/ymax;

Fx = Ex*C.q_0;
Fy = Ey*C.q_0;

ax = Fx/C.m_e;
ay = 0;

ef = sqrt(Ex^2 + Ey^2);
vth = sqrt((2 * C.kb * C.temp) / C.m_e);

tau = 0.2e-12;
mfpath= vth*tau;

x = xmax*rand(nAtom,1);
y = ymax*rand(nAtom,1);

sigma = sqrt((C.kb * C.temp) / C.m_e);
MB = makedist('Normal', 'mu', 0, 'sigma', sigma);
Vx = icdf(MB,rand(nAtom,1));
Vy = icdf(MB, rand(nAtom,1));

time = linspace(0,tStep*nStep, nStep);

color = hsv(nPlot);

scatter = 1- exp(-tStep/0.2e-12);

for i = 1:nStep

    Ux = logical(x>=xmax);
    Lx = logical(x<=0);
    
    Uy = logical(y>=ymax);
    Ly = logical(y<=0);   

    x(Ux) = 0;
    
    x(Lx) = xmax;

    Vy(Uy) = -Vy(Uy);
    y(Uy) = ymax;
    

    Vy(Ly) = -Vy(Ly);
    y(Ly) = 0;
      
    for j=1:length(x)
        if scatter > rand()
                Vx(j)= icdf(MB, rand());
                Vy(j)= icdf(MB, rand());
        end
    end
    

    yp = y;
    xp = x; 
       
    Vx = Vx + ax*tStep;
    Vy = Vy + ay*tStep;
    
    x = x + Vx*tStep;
    y = y + Vy*tStep;
    
   
    
    figure(1)
    xlabel('(m)');
    ylabel('(m)');
    title('Particle trajectory');
    axis ([0 xmax 0 ymax]);
    pause(0.0001)
    hold on;
    
    for j=1:nPlot 
        plot([xp(j)';x(j)'], [yp(j)';y(j)'], 'color', color(j,:)); 
    end 
    
    Vavg = mean(Vx.^2 + Vy.^2);
    dens(i) = sqrt(Vavg)*C.q_0*Econ;    
    mfpath = tau*Vavg; 
    tcalc = (mfpath*nAtom)/Vavg;
    

    
end

figure(2)
plot(time,dens, '.');
title('Current density over time');
xlabel('Time (s)');
ylabel('Drift Current (A/m)');

T = zeros(100,50);
NumPart = zeros(100,50);
xbins = linspace(0,xmax, 100);
ybins = linspace(0,ymax, 50);

for i = 1:nAtom
    for x_idx = 1:99
        for y_idx = 1:49
            if (x(i,1) > xbins(x_idx) && x(i,1) <= xbins(x_idx+1) && y(i,1) > ybins(y_idx) && y(i,1) <= ybins(y_idx+1))
                T(x_idx,y_idx) = T(x_idx,y_idx) + ((Vx(i,1)^2 + Vy(i,1)^2)*C.m_e)/(2*C.kb);
                NumPart(x_idx,y_idx) = NumPart(x_idx,y_idx) + 1;
            end
        end
    end
end

[x_b, y_b] = meshgrid(1:2:100,1:2:200);

Temp = T./NumPart;
Temp(isnan(Temp))=0;
figure(3);
surf(x_b, y_b, Temp);
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');
colorbar;

density = NumPart./sum(NumPart); 
density(isnan(density))=0; 
figure(4);
surf(x_b, y_b, density);
colorbar;
grid on;
title('Density Map')
xlabel('x position (nm)');
ylabel('y position (nm)');
colorbar;
hold off;