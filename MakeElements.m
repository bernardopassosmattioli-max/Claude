function [NosEle]=MakeElements(NumEleHor,NumEleVer,FirstNodePerCol)
%=================================================================%
%  Função para produzir Elementos Finitos
%  ne: matriz que fornece os nós de cada elemento
%             O----------------------O
%             |1                  2  |
%             |     ELEMENTO         |
%             |                      |
%             |4                  3  |
%             O----------------------O
%  FirstNodePerCol: Número do primeiro nó de cada coluna
%  
%=================================================================%
NosEle = zeros(NumEleHor * NumEleVer , 4);
for i=1:NumEleHor 
    NosEle([i * NumEleVer : -1:(i-1) * NumEleVer + 1] , 4) = [FirstNodePerCol(i): -1:FirstNodePerCol(i) - NumEleVer + 1].';
end;
NosEle(:,1)=NosEle(:,4)- 1;
NosEle(:,2)=NosEle(:,1)+ NumEleVer + 1;
NosEle(:,3)=NosEle(:,2)+ 1;
% THE END OF THIS FUNCTION