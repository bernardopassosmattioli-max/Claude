%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk;

K = zeros(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = zeros(2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);

for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
end

Lx = 10; Ly = 5;    % dimensoes da Estrutura 2 [m]
no = @(xp,yp) (nely+1)*round(xp/Lx*nelx) + round(yp/Ly*nely) + 1;   % no na posicao (xp,yp) [m]
F(2*no(5,0),1) = -1;      % carga em x = 5 m, y = 0
F(2*(nelx+1)*(nely+1)-1,1) = 0;
fixeddofs = [2*no(0,0), 2*no(10,0)-1, 2*no(10,0)];   % apoios em (0,0) e (10,0)

alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end