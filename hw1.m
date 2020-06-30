
[x,y] = meshgrid(-5:0.25:5,-5:0.25:5);
%a
i=cos(x)
j=sin(y)
figure
subplot(1,3,1)
q=quiver(x,y,i,j)

subplot(1,3,2)
c=curl([i,j],[x y])
sc=surf(c)
colormap(hsv)


subplot(1,3,3)
d=divergence([i,j],[x y])
sd=surf(d)
colormap(hsv)
sd.EdgeColor='none'

%b
i=cos(y)
j=sin(x)
figure
subplot(1,3,1)
q=quiver(x,y,i,j)

subplot(1,3,2)
c=curl([i,j],[x y])
sc=surf(c)
colormap(hsv)
sc.EdgeColor='none'

subplot(1,3,3)
d=divergence([i,j],[x y])
sd=surf(d)
colormap(hsv)
sd.EdgeColor='none'

%c
i=cos(x*y)
j=sin(x*y)
figure
subplot(1,3,1)
q=quiver(x,y,i,j)

subplot(1,3,2)
c=curl([i,j],[x y])
sc=surf(c)
colormap(hsv)
sc.EdgeColor='none'

subplot(1,3,3)
d=divergence([i,j],[x y])
sd=surf(d)
colormap(hsv)
sd.EdgeColor='none'

%d
i=y.^2
j=x.^2
figure
subplot(1,3,1)
q=quiver(x,y,i,j)

subplot(1,3,2)
c=curl([i,j],[x y])
sc=surf(c)
colormap(hsv)
sc.EdgeColor='none'

subplot(1,3,3)
d=divergence([i,j],[x y])
sd=surf(d)
colormap(hsv)
sd.EdgeColor='none'