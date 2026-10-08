function [TVM, Stress] = VonMises(nelx,nely,U,E,nu,a,b)

% ================================================================
% CALCULO DA TENSAO DE VON MISES NO CENTRO DE CADA ELEMENTO
%
% nelx  = numero de elementos na horizontal
% nely  = numero de elementos na vertical
% U     = vetor global de deslocamentos
% E     = modulo de Young
% nu    = coeficiente de Poisson
% a     = largura do elemento
% b     = altura do elemento
%
% TVM(el) = tensao de von Mises do elemento
%
% ================================================================


%% MATRIZ CONSTITUTIVA C
% Estado plano de tensao

C = E/(1-nu^2) * ...
    [ 1    nu   0;
      nu   1    0;
      0    0   (1-nu)/2 ];


%% MATRIZ V PARA VON MISES

V = [ 1    -1/2   0;
     -1/2   1     0;
      0     0     3 ];


%% MATRIZ B NO CENTRO DO ELEMENTO
%
% Elemento quadrilateral de 4 nos:
%
%      4 -------- 3
%      |          |
%      |          |
%      |          |
%      1 -------- 2
%

% Derivadas das funcoes de forma no centro
% xi = 0
% eta = 0

dNdx = [ -1/(2*a), ...
          1/(2*a), ...
          1/(2*a), ...
         -1/(2*a) ];

dNdy = [ -1/(2*b), ...
         -1/(2*b), ...
          1/(2*b), ...
          1/(2*b) ];


% Matriz B

B = zeros(3,8);

for i = 1:4

    B(1,2*i-1) = dNdx(i);

    B(2,2*i)   = dNdy(i);

    B(3,2*i-1) = dNdy(i);
    B(3,2*i)   = dNdx(i);

end


%% INICIALIZACAO

NEle = nelx*nely;

TVM = zeros(NEle,1);

Stress = zeros(NEle,3);


%% PERCORRE TODOS OS ELEMENTOS

el = 0;

for elx = 1:nelx

    for ely = 1:nely

        el = el + 1;


        %% IDENTIFICACAO DOS NOS DO ELEMENTO

        n1 = (nely+1)*(elx-1) + ely;

        n2 = (nely+1)*elx + ely;


        %% GRAUS DE LIBERDADE DO ELEMENTO

        edof = [2*n1-1;
                2*n1;
                2*n2-1;
                2*n2;
                2*n2+1;
                2*n2+2;
                2*n1+1;
                2*n1+2];


        %% VETOR DE DESLOCAMENTOS DO ELEMENTO

        ue = U(edof);


        %% DEFORMACOES

        strain = B*ue;


        %% TENSOES

        sigma = C*strain;


        %% SALVA AS TENSOES

        Stress(el,:) = sigma';


        %% TENSAO DE VON MISES

        TVM(el) = sqrt(sigma' * V * sigma);

    end

end

end