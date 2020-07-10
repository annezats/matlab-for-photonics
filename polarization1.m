%% Assignment Instructions
% Write code in Matlab that calculates the polarization vector as a function
% of position for linearly, circularly, and elliptically polarized light
% propagating in the z direction plot and its projection onto the x-y axis.
% Visualize this data to reproduce the plots on slide 17

%% Constants
t = linspace(0,40); %independent variable for the wave function
 
E1=5; %Amplitude
E2=5;

k1= 0.1*pi; %Frequency
k2= 0.1*pi;

%% Linearly
omega1=4*pi/2;
omega2=4*pi/2;
x1 = E1.*cos(k1*t + omega1);
y1 = E2.*cos(k2*t + omega2);
figure
plot3(t,x1,y1,'LineWidth',1.5);
title('Linear Polarization')
xlabel('Z'); 
ylabel('X');
zlabel('Y');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);

hold on
plot3(t*0,x1,y1, 'LineWidth',1.5);

axis equal
xl = xlim();

hold on 
line(xl, [0,0], [0,0], 'LineWidth', 2, 'Color', 'k');



%% Circularly
omega3=4*pi/2; %Phase
omega4=1*pi/2;
x2 = E1.*cos(k1*t + omega3);
y2 = E2.*cos(k2*t + omega4);
figure
plot3(t,x2,y2,'LineWidth',1.5);
title('Circular Polarization')
xlabel('Z'); 
ylabel('X');
zlabel('Y');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);
hold on
plot3(t*0,x2,y2, 'LineWidth',1.5);

axis equal
xl = xlim();

hold on 
line(xl, [0,0], [0,0], 'LineWidth', 2, 'Color', 'k');


%% Elliptically

omega5=4*pi/2;
omega6=pi/6;
x3 = E1.*cos(k1*t + omega5);
y3 = E2.*cos(k2*t + omega6);
figure
plot3(t,x3,y3,'LineWidth',1.5);
title('Elliptical Polarization')
xlabel('Z'); 
ylabel('X');
zlabel('Y');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);

hold on
plot3(t*0,x3,y3, 'LineWidth',1.5);

axis equal
xl = xlim();

hold on 
line(xl, [0,0], [0,0], 'LineWidth', 2, 'Color', 'k');


