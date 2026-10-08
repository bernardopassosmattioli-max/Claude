clear; clc;

nelx = 10;
nely = 24;
rho(1:nely, 1:nelx) = 1;

a = 1; b = 1;
M = [nely + 1, nelx + 1];
[x, y] = meshgrid(a * [0 : nelx], b * [0 : nely]);
FirstNodePerCol = find(y(:) == max(y(:)));
NosX = -x(:); NosY = y(:);
NosEle = MakeElements(nelx, nely, FirstNodePerCol);
NEle = (nelx) * (nely);
Node = [NosX NosY];
NNos = (nelx + 1) * (nely + 1);

[U] = FE(nelx, nely);

[KE] = lk;

E = 1;

nu = 0.3;

[TVM, Stress] = VonMises(nelx,nely,U,E,nu,a,b);

figure(1); colormap(gray); imagesc(-rho); axis equal; axis tight; axis off;

for i = 1: 0.5 * length(U)
    Ux(i,1) = U(2 * i - 1);
    Uy( i ,1) = U(2 * i );
end
figure(2); surf(reshape(Ux, [nely+1, nelx+1])); set(gca, 'XDir','normal');
figure(3); surf(reshape(Uy, [nely+1, nelx+1])); set(gca, 'YDir','reverse');

TVM_Matrix = reshape(TVM,[nely,nelx]);

figure

imagesc(TVM_Matrix)

axis equal
axis tight

colorbar

title('Tensao de Von Mises')

xlabel('Elementos na direcao X')

ylabel('Elementos na direcao Y')
