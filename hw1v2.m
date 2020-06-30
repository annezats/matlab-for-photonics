[x,y] = meshgrid(-5:0.25:5,-5:0.25:5);
%a
i=cos(x)
j=sin(y)
figure
subplot(1,3,1)
quiver(x,y,i,j)
subplot(1,3,2)
cav = curl(x,y,i,j)
pcolor(x,y,cav);
subplot(1,3,3)
d=divergence(x,y,i,j)
pcolor(x,y,d)
%b
i=cos(y)
j=sin(x)
figure
subplot(1,3,1)
quiver(x,y,i,j)
subplot(1,3,2)
cav = curl(x,y,i,j)
pcolor(x,y,cav);
subplot(1,3,3)
d=divergence(x,y,i,j)
pcolor(x,y,d)
%c
i=cos(x*y).
j=sin(x*y)
figure
subplot(1,3,1)
quiver(x,y,i,j)
subplot(1,3,2)
cav = curl(x,y,i,j)
pcolor(x,y,cav);
subplot(1,3,3)
d=divergence(x,y,i,j)
pcolor(x,y,d)
%d
i=y.^2
j=x.^2
figure
subplot(1,3,1)
quiver(x,y,i,j)
subplot(1,3,2)
cav = curl(x,y,i,j)
pcolor(x,y,cav);
subplot(1,3,3)
d=divergence(x,y,i,j)
pcolor(x,y,d)