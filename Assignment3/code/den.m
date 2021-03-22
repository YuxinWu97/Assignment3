function [J_dens, calcurrent] = den(Lb, Wb)
global C
C.temp = 300;                       % Initial temperature 
C.kb = 1.3806504e-23;               % Boltzmann constant
C.m_0 = 9.10938215e-31;             % Electron mass
C.m_e = 0.26*C.m_0;                 % Effective mass 
C.q_0 = 1.60217662e-19;               % Charge of electron




si = 1;          %Initialize resistivity inside the bottleneck
so = 10e-2;      %Initialize resistivity outside the bottleneck


[Current, V_solution, Ex, Ey] = cur(200,100, Lb, Wb, si, so);

calcurrent = Current; 

[X,Y] = meshgrid(1:100,1:200);

Econ = 1e15;

Econ = Econ/0.001;

xmax = 200e-9;
ymax = 100e-9;

x1 = xmax/2 - Lb;
x2 = xmax/2 + Lb;

y1 = ymax/2 - Wb;
y2 = ymax/2 + Wb;

tStep = .01e-12;
nStep = 50;
nAtom = 3000;
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
    while (y(i)<= y1 || y(i)>= y2) && (x(i)>= x1 && x(i) <= x2)
        x(i) = xmax*rand(); 
        y(i) = ymax*rand(); 
    end
end


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
    
     if specular == 1 
         for j=1:nAtom
            if (y(j)<= y1 || y(j) >= y2) && (x(j)>= x1 && x(j) <= x2)
                Vx(j) = - Vx(j);
                x(j) = xp(j);
                y(j) = yp(j);
            end
            if (y(j) <= y1 && y(j) >= y2) && (x(j) >= x1 && x(j) <= x2)
                Vy(j) = - Vy(j);
                x(j) = xp(j);
                y(j) = yp(j);
            end
         end
     else 
            for j=1:nAtom
                while (y(j)<= y1 || y(j)>= y2) && (x(j)>= x1 && x(j) <= x2)
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
        for x_idx = 1:200
            for y_idx = 1:100
                if (x(j,1) > xbins(x_idx) && x(j,1) <= xbins(x_idx+1) && y(j,1) > ybins(y_idx) && y(j,1) <= ybins(y_idx+1))
                        Vx(j,1) = Vx(j,1) + (ax(x_idx, y_idx)/1e-9)*tStep;
                        Vy(j,1) = Vy(j,1) + (ay(x_idx, y_idx)/1e-9)*tStep;
                       
                end
            end
        end
    end
      
    x = x + Vx*tStep;
    y = y + Vy*tStep;

    Vavg = mean(Vx.^2 + Vy.^2);

    J_dens_vect(i) = sqrt(Vavg)*C.q_0*Econ;

    
end

J_dens = mean(J_dens_vect);


end
