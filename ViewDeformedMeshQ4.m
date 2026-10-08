NumGL = 2;
CoordDef=zeros(NNos,2); CoordOrig = CoordDef;
Uf=zeros(NNos,NumGL);

c=1;
b=1;
for i=1:NNos
    for j=1:NumGL
        Uf(i,j)=U(b);
        b=b+1;
    end
end

% Fator de amplificacao para visualizacao
fator = 0.1;

for ii=1:NNos
  	CoordDef(ii,1)=-(Node(ii,1) + fator*Uf(ii,1));
  	CoordDef(ii,2)=Node(ii,2) + fator*Uf(ii,2);
end
  
Mesh = ViewMeshQ4(NNos, NEle, NosEle);
figure(4); gplot(Mesh,CoordDef,'-k'); axis equal;

% Plotar malha original em vermelho
for ii=1:NNos
  	CoordOrig(ii,1)=-Node(ii,1);
  	CoordOrig(ii,2)=Node(ii,2);
end
hold on; gplot(Mesh,CoordOrig,'-r'); axis equal;
