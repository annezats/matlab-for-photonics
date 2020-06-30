

[x,y] = meshgrid(-5:0.5:5,-5:0.5:5);
I=[cos(x)  cos(y)  cos(x*y)  y.^2]
J=[sin(y)  sin(x)  sin(x*y)  x.^2]
for a=1:4
    i=I(a)
    j=J(a)
    figure
    subplot(1,3,1)
    q=quiver(x,y,i,j)

    hold on
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
end
%fi=diff(i)
%fj=diff(j)
%diverg= dot([fi,fj],[i,j]) %this is how divergence works


