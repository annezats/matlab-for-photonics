[x,y] = meshgrid(-5:1:5,-5:1:5);
%a
i=cos(x)
j=sin(y)
figure
quiver(x,y,i,j) 
%b
i=cos(y)
j=sin(x)
figure 
quiver(x,y,i,j)
%c
i=cos(x*y)
j=sin(x*y)
figure 
quiver(x,y,i,j)
%d
i=y.^2
j=x.^2
figure 
quiver(x,y,i,j)
