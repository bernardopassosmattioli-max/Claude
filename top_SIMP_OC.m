%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
clear; clc;
% DADOS DE ENTRADA
estrutura = 1;   % 1 = Estrutura 1, 2 = Estrutura 2, 3 = Exercicio 2
if estrutura==1, nelx=40; nely=96;
elseif estrutura==2, nelx=80; nely=40;
else nelx=50; nely=30;
end
h = 10/nelx;     % tamanho do elemento [m]
volfrac = 0.3;
penal = 3;
rmin = 1.5;
% REGIAO PASSIVA (recorte vazio do Exercicio 2)
passive = zeros(nely,nelx);
if estrutura==3
    for elx = 1:nelx
        for ely = 1:nely
            xc = (elx-0.5)*h; yc = (ely-0.5)*h;
            if xc>3 && xc<7 && yc<3, passive(ely,elx) = 1; end
        end
    end
end
% INITIALIZE
x(1:nely,1:nelx) = volfrac;
x(passive==1) = 0.001;
loop = 0;
change = 1.;
historico = [];
% START ITERATION
while change > 0.01
    loop = loop + 1;
    xold = x;
    % FE-ANALYSIS
    if estrutura==1, U = FE(nelx,nely,x,penal);
    elseif estrutura==2, U = FEb(nelx,nely,x,penal);
    else U = FE2(nelx,nely,x,penal);
    end
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk;
    c = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx   +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
            c = c + x(ely,elx)^penal*Ue'*KE*Ue;
            dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
        end
    end
    historico(loop) = c;
    % FILTERING OF SENSITIVITIES
    [dc]   = check(nelx,nely,rmin,x,dc);
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [x]    = OC(nelx,nely,x,volfrac,dc,passive);
    % PRINT RESULTS
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
    % PLOT DENSITIES
    colormap(gray); imagesc(-flipud(x)); axis equal; axis tight; axis off;pause(1e-3);
end
% CURVA DE CONVERGENCIA
figure; plot(1:loop,historico,'-o'); grid on;
xlabel('Iteracao'); ylabel('Funcao objetivo c'); title('Convergencia');
% TENSAO DE VON MISES
[TVM] = VonMises(nelx,nely,U,1,0.3,h,h);
TVM = TVM.*x(:).^penal;
figure; imagesc(flipud(reshape(TVM,nely,nelx))); axis equal; axis tight; colorbar;
title('Tensao de von Mises');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
