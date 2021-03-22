%Part2
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

si= 1;
so= 10e-2;


[Current, V_solution, Ex, Ey] = cur (200,100, 40, 40, si, so);


[X,Y] = meshgrid(1:100,1:200);

figure();
surf(X,Y,V_solution);
grid;
title('Plot of V(x,y)');
xlabel('y (nm)');
ylabel('x (nm)');
zlabel('Voltage (A)');

figure();
quiver(Ex', Ey');
axis([0 200 0 100]);
grid;
title('Electric field vector');
xlabel('x (nm)');
ylabel('y (nm)');




xmax = 200e-9;
ymax = 100e-9;

tStep = .01e-12;

nStep = 1000;
nAtom = 1000;
nPlot = 50;

specular = 1;

T_avg = zeros(nStep,1);
time = linspace(0,tStep*nStep, nStep);


Fx = Ex*C.q_0;
Fy = Ey*C.q_0;

ax = Fx/C.m_e;
ay = Fy/C.m_e;



ef = mean(sqrt(Ex.^2 + Ey.^2));

v_th = sqrt((2 * C.kb * C.temp) / C.m_e);

sigma = sqrt((C.kb * C.temp) / C.m_e);
MB = makedist('Normal', 'mu', 0, 'sigma', sigma);
Vx = icdf(MB,rand(nAtom,1));
Vy = icdf(MB, rand(nAtom,1));

tau = 0.2e-12;
mfpath= v_th*tau;

x = xmax*rand(nAtom,1);
y = ymax*rand(nAtom,1);


for i=1:nAtom
    while (y(i)<= 40e-9 || y(i)>= 60e-9) && (x(i)>= 80e-9 && x(i) <= 120e-9)
        x(i) = xmax*rand(); 
        y(i) = ymax*rand(); 
    end
end 


color = hsv(nPlot);

scatter = 1- exp(-tStep/0.2e-12);


figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
hold on;
rectangle('position', [80e-9 0 40e-9 40e-9]);
rectangle('position', [80e-9 60e-9 40e-9 40e-9]);
hold on;
xlim(axes1,[0 xmax]);
ylim(axes1,[0 ymax]);
xlabel('nm');
ylabel('nm');

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
    
     if specular == 1 
         for j=1:nAtom
            if (y(j)<= 40e-9 || y(j) >= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                Vx(j) = - Vx(j);
                x(j) = xp(j);
                y(j) = yp(j);
            end
            if (y(j) <= 40e-9 && y(j) >= 60e-9) && (x(j) >= 80e-9 && x(j) <= 120e-9)
                Vy(j) = - Vy(j);
                x(j) = xp(j);
                y(j) = yp(j);
            end
         end
     else
            for j=1:nAtom
                while (y(j)<= 40e-9 || y(j)>= 60e-9) && (x(j)>= 80e-9 && x(j) <= 120e-9)
                    x(j) = xmax*rand(); 
                    y(j) = ymax*rand();
                    x(j) = xp(j);
                    y(j) = yp(j);
                end
            end    
     end
     
     
    yp = y;
    xp = x; 
    
    xbins = (0:1:200)*1e-9;
    ybins = (0:1:200)*1e-9;
    
    for j = 1:nAtom
        for xidx = 1:200
            for yidx = 1:100
                if (x(j,1) > xbins(xidx) && x(j,1) <= xbins(xidx+1) && y(j,1) > ybins(yidx) && y(j,1) <= ybins(yidx+1))
                        Vx(j,1) = Vx(j,1) + (ax(xidx, yidx)/1e-9)*tStep;
                        Vy(j,1) = Vy(j,1) + (ay(xidx, yidx)/1e-9)*tStep;
                       
                end
            end
        end
    end
    
    x = x + Vx*tStep;
    y = y + Vy*tStep;

    xlabel('(m)');
    ylabel('(m)');
    title('Particle trajectory');
    axis ([0 xmax 0 ymax]);
    hold on;
       
    for j=1:nPlot 
        
        plot([xp(j)';x(j)'], [yp(j)';y(j)'], 'color', color(j,:)); 
            
    end 
    
    Vavg = mean(Vx.^2 + Vy.^2); 
    
    J_dens(i) = sqrt(Vavg)*C.q_0*Econ;
    
    mfpath = tau*Vavg;
    tcalc = (mfpath*nAtom)/Vavg;

    
end

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