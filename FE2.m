%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE2(nelx,nely,x,penal)
[KE] = lk;
% Inicialização
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

% DEFINE LOADS AND SUPPORTS
%HALF MBB-BEAM
% F(2,1) = -1;
% fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
%EXERCICIO 2
Lx = 10; Ly = 6;    % dimensoes do dominio do Exercicio 2 [m]
no = @(xp,yp) (nely+1)*round(xp/Lx*nelx) + round(yp/Ly*nely) + 1;   % no na posicao (xp,yp) [m]
F(2*no(10,0),1) = -1;     % carga em x = 10 m, y = 0
F(2*(nelx+1)*(nely+1)-1,1) = 0;
fixeddofs = [1:2*no(0,2)];   % engaste em x = 0, de y = 0 ate y = 2 m
%
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end