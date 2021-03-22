function [Current, V_solution, Ex, Ey] = cur(L,W, Lb, Wb, si, so)
G = sparse(W*L,W*L);
F = sparse(W*L,1);
sigma = zeros(L, W);

dx = 1;
dy = 1;
V0 = 5;

for i=1:L
    for j=1:W

        in_x = logical( i >= (L- Lb)/2 && i <= (L + Lb)/2);
        in_y = logical( j <= Wb | j >= (W-Wb));

        if in_x && in_y
            sigma(i,j) = so;
        else 
            sigma(i,j) = si;
        end 
        
        
    end
end


[X,Y] = meshgrid(1:W,1:L);

for i = 1:L
    for j = 1:W
        n = j + (i - 1) * W;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            F(n,1) = 1;
        elseif i == L
            G(n, :) = 0;
            G(n, n) = 1; %V0 = 1
        elseif j == 1
            nxm = j + (i - 2) * W;
            nxp = j + (i) * W;
            nyp = j + 1 + (i - 1) * W;

            rxm = (sigma(i, j) + sigma(i - 1, j)) / 2.0;
            rxp = (sigma(i, j) + sigma(i + 1, j)) / 2.0;
            ryp = (sigma(i, j) + sigma(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  W
            nxm = j + (i - 2) * W;
            nxp = j + (i) * W;
            nym = j - 1 + (i - 1) * W;

            rxm = (sigma(i, j) + sigma(i - 1, j)) / 2.0;
            rxp = (sigma(i, j) + sigma(i + 1, j)) / 2.0;
            rym = (sigma(i, j) + sigma(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*W;
            nxp = j + (i)*W;
            nym = j-1 + (i-1)*W;
            nyp = j+1 + (i-1)*W;

            rxm = (sigma(i,j) + sigma(i-1,j))/2.0;
            rxp = (sigma(i,j) + sigma(i+1,j))/2.0;
            rym = (sigma(i,j) + sigma(i,j-1))/2.0;
            ryp = (sigma(i,j) + sigma(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end


V = G\F; 

V_solution = zeros(L,W);

for i=1:L
    for j=1:W
        n = j + (i-1)*W;
        
        V_solution(i,j) = V(n);
    end
end
    

[X,Y] = meshgrid(1:W,1:L);

for i = 1:L
    for j = 1:W
        if i == 1
            Ex(i, j) = (V_solution(i + 1, j) - V_solution(i, j));
        elseif i == L
            Ex(i, j) = (V_solution(i, j) - V_solution(i - 1, j));
        else
            Ex(i, j) = (V_solution(i + 1, j) - V_solution(i - 1, j)) * 0.5;
        end
        
        if j == 1
            Ey(i, j) = (V_solution(i, j + 1) - V_solution(i, j));
        elseif j == W
            Ey(i, j) = (V_solution(i, j) - V_solution(i, j - 1));
        else
            Ey(i, j) = (V_solution(i, j + 1) - V_solution(i, j - 1)) * 0.5;
        end
    end
end


Ex = -Ex;
Ey = -Ey; 

Fx = sigma .* Ex;
Fy = sigma .* Ey;

C0 = sum(Fx(1, :));
Cnx = sum(Fx(L, :));
Current = (C0 + Cnx) * 0.5;