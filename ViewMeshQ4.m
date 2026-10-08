% Gera matriz M para plotagem
% de malhas quadrilaterais de 4 nos

function [M] = ViewMeshQ4(NNos, NEle, NosEle)

M = zeros(NNos,NNos);
for e = 1:NEle
   noi = NosEle(e,1);               
   noj = NosEle(e,2);               
   nok = NosEle(e,3);
   nol = NosEle(e,4);
   
   M(noi,noi)=1;
   M(noi,noj)=1;
   M(noi,nok)=1;
   M(noi,nol)=1;
   
   M(noj,noi)=1;
   M(noj,noj)=1;
   M(noj,nok)=1;
   M(noj,nol)=1;
   
   M(nok,noi)=1;
   M(nok,noj)=1;
   M(nok,nok)=1;
   M(nok,nol)=1;
   
   M(nol,noi)=1;
   M(nol,noj)=1;
   M(nol,nok)=1;
   M(nol,nol)=1;
end

