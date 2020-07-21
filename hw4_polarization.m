%calculate the polarization vector as a function of position 
%for linearly, circularly, and elliptically polarized light
%propogating in the z direction, and its projection onto the x-y axis
clear

%some constants
k= pi/2 %wave vector (frequency)
%our direction of propogation
z=0:0.1:2*pi

%% Linear Polarization
E1= 1 %amplitude
E2= 1
o1=0 %(omega/angular freq*time=phase (changes w time))
o2=0 %the two curves being in phase is what makes it linear
x=E1*sin(k*z-o1);
y=E2*sin(k*z-o2);

figure 
hold on
%----total polarization
plot3(z,x,y,'LineWidth',2)
%----x-y projection
plot3(z*0+2*pi,x,y,'LineWidth',2)
%-----x and y components
plot3(z,x,y*0,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x,y*0,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
plot3(z,x*0,y,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x*0,y,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
%----making a line through the origin
plot(z,z*0,'LineWidth',2)
%---making it pretty
grid on
xlim([0 2*pi]);
ylim([-1*E1 E1]);
zlim([-1*E2 E2]);
title('linear polarization')
xlabel('z'); 
ylabel('x');
zlabel('y');
view(-60,40);
hold off

%% Circular Polarization
E3= 1 %amplitude
E4= 1
o3=0 %(omega/angular freq*time=phase (changes w time))
o4=pi/2 %the two curves being 90 degrees out of phase is what makes it circular
x=E3*sin(k*z-o3);
y=E4*sin(k*z-o4);

figure 
hold on
%----total polarization
plot3(z,x,y,'LineWidth',2)
%----x-y projection
plot3(z*0+2*pi,x,y,'LineWidth',2)
%-----x and y components
plot3(z,x,y*0,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x,y*0,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
plot3(z,x*0,y,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x*0,y,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
%----making a line through the origin
plot(z,z*0,'LineWidth',2)
%---making it pretty
grid on
xlim([0 2*pi]);
ylim([-1*E3 E3]);
zlim([-1*E4 E4]);
title('circular polarization')
xlabel('z'); 
ylabel('x');
zlabel('y');
view(-60,40);
hold off

%% Elliptical Polarization
E5= 1 %amplitude
E6= 5
o5=0 %(omega/angular freq*time=phase (changes w time))
o6=pi/3 
x=E5*sin(k*z-o5);
y=E6*sin(k*z-o6);

figure 
hold on
%----total polarization
plot3(z,x,y,'LineWidth',2)
%----x-y projection
plot3(z*0+2*pi,x,y,'LineWidth',2)
%-----x and y components
plot3(z,x,y*0,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x,y*0,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
plot3(z,x*0,y,'LineWidth',2)
quiver3(z,0*z,0*z,0*z,x*0,y,'AutoScale','off','LineWidth',0.5,'ShowArrowHead', 'off')
%----making a line through the origin
plot(z,z*0,'LineWidth',2)
%---making it pretty
grid on
xlim([0 2*pi]);
ylim([-1*E5 E5]);
zlim([-1*E6 E6]);
title('elliptical polarization')
xlabel('z'); 
ylabel('x');
zlabel('y');
view(-60,40);
hold off

%Brought to you by Anne Zats, Natascha Krishnanand, and Anton Kyrylenko