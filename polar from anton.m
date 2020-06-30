%Polarization 

t=0:.001:40; %independent variable for the polarization
s=linspace(0,40); %independent variable for the wave function

%constants for equations 
E1=1; %Amplitude
E2=1;
omega1=4*pi/2; %Phase
omega2=1*pi/2;
k1=0.4*pi; %Frequency
k2=0.4*pi;



%ciruclar polarization
x1 = E1.*cos(k1*t + omega1);
y1 = E2.*cos(k2*t + omega2);
figure
plot3(t,x1,y1,'LineWidth',1.5);
title('circular polarization')
xlabel('wave direction'); 
ylabel('H');
zlabel('E');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);



%elliptical polarization
omega3=4*pi/2;
omega4=pi/6;
x2 = E1.*cos(k1*t + omega3);
y2 = E2.*cos(k2*t + omega4);
figure
plot3(t,x2,y2, 'LineWidth',1.5);
title('elliptical polarization')
xlabel('wave direction'); 
ylabel('H');
zlabel('E');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);



%linear polarization
omega5=4*pi/2;
omega6=4*pi/2;
x3 = E1.*cos(k1*t + omega5);
y3 = E2.*cos(k2*t + omega6);
figure
plot3(t,x3,y3, 'LineWidth',1.5);
hold on
d = ones(size(s));
c = E1*cos(k1*s+omega5*d);
patch([s fliplr(s)], [c zeros(size(c))], [zeros(size(s)) zeros(size(s))], 'c')
patch([s fliplr(s)], [zeros(size(s)) zeros(size(s))], [c zeros(size(c))], 'y')
hold off
title('linear polarization')
xlabel('wave direction'); 
ylabel('H');
zlabel('E');
xlim([-5 50]);
ylim([-5 5]);
zlim([-5 5]);
grid on
set(gca,'FontSize',14);