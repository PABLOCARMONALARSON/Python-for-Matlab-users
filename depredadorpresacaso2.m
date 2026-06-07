function depredadorpresacaso2
model = createpde(2);
L=500;
gdm = [3;4;0;L;L;0;0;0;L;L];
g = decsg(gdm,'R1',('R1')');
geometryFromEdges(model,g);
gm=generateMesh(model,'Hmax', 5);
pdemesh(gm)
pause
%constantes del modelo adimensional
D1=1; 
D2=10; 
r=1; 
K=1; 
b=2.5; 
e=2.5; 
m=5;
dd=0.6; 

d=[1;1];
c=[D1;D2];

f = @(location,state) [r*state.u(1,:).*(1-state.u(1,:)./K)-b*(state.u(1,:).*state.u(2,:))./(1+e*state.u(1,:));...
    m*(state.u(1,:).*state.u(2,:))./(1+e*state.u(1,:))-dd*state.u(2,:)];
                   
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',0,'f',f);

endTime = 150;
tlist = 0:0.5:endTime;
setInitialConditions(model,@u0);
g=[0;0];
applyBoundaryCondition(model,'neumann','Edge',1:4,'g',g,'q',g);
%Resolucion del problema
disp('Resolucion')
R = solvepde(model,tlist);
disp('Acabada')
u = R.NodalSolution;
msh = R.Mesh;
figure;
for i = 1:length(tlist)
    if(mod(i,20)==0)
    subplot(2,1,1)
    pdeplot(msh,XYData=u(:,1,i),ZData=u(:,1,i));
    %pdeplot(model, 'XYData', u(:,1,i), 'Contour', 'off');
    title(['Presas en t = ' num2str(tlist(i))]);
    axis equal;
    colormap jet;
    %caxis([0 0.7]); % Fijar escala de color
    drawnow;
    subplot(2,1,2)
    pdeplot(msh,XYData=u(:,2,i),ZData=u(:,2,i));
   % pdeplot(model, 'XYData', u(:,2,i), 'Contour', 'off');
    title(['Depredadores en t = ' num2str(tlist(i))]);
    axis equal;
    colormap jet;
    %caxis([0 1.3]); % Fijar escala de color
    drawnow;
    end
    
end
 

function ui=u0(location,state)
yb=location.y;
xb=location.x;  
ui =[6./35-2.e-7*(xb-180).*(xb-720)-6.e-7*(yb-90).*(yb-210);...
   116./245-3.e-5*(xb-450)-6e-5*(yb-135)];
end
end