%%%%%%%%%% FE-ANALYSIS - EXERCICIO 2 (dominio 10 x 6 com recorte) %%%%%%%%%%
function [U]=FE2(nelx,nely,x,penal)
[KE] = lk;
h = 10/nelx;   % tamanho do elemento [m]

K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
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

% Carga: canto inferior direito, para baixo
F(2*((nely+1)*nelx + 1),1) = -1;
% Engaste: borda esquerda, de y = 0 ate y = 2
fixeddofs = [1:2*(round(2/h)+1)];

alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end
