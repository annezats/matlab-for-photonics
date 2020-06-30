%getting a pretty wave
%calculate the polarization vector as a function of position 
%for linearly, circularly, and elliptically polarized light
%propogating in the z direction, and its projection onto the x-y axis
%visualize this
clear

%some constants
E0= 1 %magnitude
E1= 1
k= pi/2 %wave vector (frequency)
o1=0 %(omega/angular freq*time=phase (changes w time))
o2=0
%our direction of propogation
z=0:0.1:2*pi

%linear polarization
x=E0*cos(k*z-o1);
E1=0
y=E1*sin(k*z-o2);

figure 
plot3(z,x,y)
hold on
plot3(z*0+2*pi,x,y)% x-y projection
%----making a line through the origin
plot(z,z*0)
%---making it go to the point
quiver3(z,0*z,0*z,0*z,x,y,'AutoScale','off')
%Plot3AxisAtOrigin(x,y,z)
%---making it pretty
grid on
xlim([0 2*pi]);
ylim([-1 1]);
zlim([-1 1]);
title('linear polarization')
xlabel('z'); 
ylabel('x');
zlabel('y');
hold off

%circular polarization
x=E0*cos(k*z-o1);
E1=E0
y=E1*sin(k*z-o2);

figure 
plot3(z,x,y)
hold on
plot3(z,z*0,y)% x wave
plot3(z,x,z*0)% y wave
plot3(z*0+2*pi,x,y)% x-y projection
hold on
%----making a line through the origin
plot(z,z*0)
%---making it go to the point
quiver3(z,0*z,0*z,0*z,x,y,'AutoScale','off')
%Plot3AxisAtOrigin(x,y,z)
grid on
xlim([0 2*pi]);
ylim([-1 1]);
zlim([-1 1]);
title('circular polarization')
xlabel('z'); 
ylabel('x');
zlabel('y');
hold off

%linear polarization
x=E0*cos(k*z-o1);
E1=2
%do I want to change the omega to get the elliptical?
y=E1*sin(k*z-o2);

figure 
plot3(z,x,y)
hold on
plot3(z,x*0,y)% x wave
plot3(z,x,y*0)% y wave
plot3(z*0+2*pi,x,y)% x-y projection
hold on
%----making a line through the origin
plot(z,z*0)
%---making it go to the point
quiver3(z,0*z,0*z,0*z,x,y,'AutoScale','off')
%Plot3AxisAtOrigin(x,y,z)
grid on
xlim([0 2*pi]);
ylim([-2 2]);
zlim([-2 2]);
title('elliptical polarization')
xlabel('t'); 
ylabel('x');
zlabel('y');
hold off


%things to do:
%rotate the axes
%--rotating the data is easier
%make a grid through zero in the x and y directions
%--copied plot3axisatorigin
%overlay the x,y view of it
%--graph the y z and set the x to 0+ offset
%fill in the middle?
%--used quiver
%make the formulas actually the right ones
%--I am so confused w the variables, nvm i kinda figured it out? we r
%ignoring time
%make 3 separate plots for linear circular n elliptical
%--done