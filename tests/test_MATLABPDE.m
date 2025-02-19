g = @scatterg;
k = 60;
c = 1;
a = -k^2;
f = 0;

numberOfPDE = 1;
pdem = createpde(numberOfPDE);

geometryFromEdges(pdem,g);

specifyCoefficients(pdem,'m',0,'d',0,'c',c,'a',a,'f',f);

figure;
pdegplot(pdem,'EdgeLabels','on');
axis equal
title 'Geometry With Edge Labels Displayed';
ylim([0,1])


bOuter = applyBoundaryCondition(pdem,'Edge',(5:8),'g',0,'q',-60i);
innerBCFunc = @(loc,state)-exp(-1i*k*loc.x);
bInner = applyBoundaryCondition(pdem,'Edge',(1:4),'u',innerBCFunc);

generateMesh(pdem,'Hmax',0.02);
figure
pdemesh(pdem);
axis equal

result = solvepde(pdem);
u = result.NodalSolution;

figure
pdeplot(pdem,'xydata',real(u),'zdata',real(u),'mesh','off');
colormap(cool)

figure
m = 10;
h = newplot;
hf = h.Parent;
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(u));
for j = 1:m
    uu = real(exp(-j*2*pi/m*sqrt(-1))*u);
    pdeplot(pdem,'xydata',uu,'colorbar','off','mesh','off');
    caxis([-maxu maxu]);
    axis tight
    ax = gca;
    ax.DataAspectRatio = [1 1 1];
    axis off
    M(j) = getframe(hf);
end
movie(hf,M,2);