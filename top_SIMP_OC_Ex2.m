%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
clear; clc;
% DADOS DE ENTRADA
nelx = 50;
nely = 30;
penal = 3;
rmin = 1.5;
h = 10/nelx;   % tamanho do elemento [m]
% REGIAO PASSIVA: recorte vazio de x = 3 a 7 m e y = 0 a 3 m
passive = zeros(nely,nelx);
for elx = 1:nelx
    for ely = 1:nely
        xc = (elx-0.5)*h; yc = (ely-0.5)*h;   % centro do elemento
        if xc > 3 && xc < 7 && yc < 3
            passive(ely,elx) = 1;
        end
    end
end
for volfrac = [0.2 0.35 0.6]
% INITIALIZE
x = volfrac*ones(nely,nelx);
x(passive==1) = 0.001;
dc = zeros(nely,nelx);
loop = 0;
change = 1.;
historico = [];
% START ITERATION
while change > 0.01
    loop = loop + 1;
    xold = x;
    % FE-ANALYSIS
    [U] = FE2(nelx,nely,x,penal);
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
    dc(passive==1) = 0;
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [x]    = OC(nelx,nely,x,volfrac,dc);
    % PRINT RESULTS
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
    % PLOT DENSITIES
    figure(1); colormap(gray); imagesc(-flipud(x)); axis equal; axis tight; axis off;pause(1e-3);
end
% TOPOLOGIA FINAL
figure; colormap(gray); imagesc(-flipud(x)); axis equal; axis tight; axis off;
title(['Topologia - ' num2str(100*volfrac) '% de material']);
% CURVA DE CONVERGENCIA
figure; plot(1:loop,historico,'-o'); grid on;
xlabel('Iteracao'); ylabel('Funcao objetivo c'); title(['Convergencia - ' num2str(100*volfrac) '% de material']);
% TENSAO DE VON MISES
[TVM] = VonMises(nelx,nely,U,1,0.3,h,h);
TVM = TVM.*x(:).^penal;
figure; imagesc(flipud(reshape(TVM,nely,nelx))); axis equal; axis tight; colorbar;
title(['Tensao de von Mises - ' num2str(100*volfrac) '% de material']);
end


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
