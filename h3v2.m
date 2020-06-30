clear
% angle of incidence (arbitrary constant)
ti=pi/4  
%-------------------------
ni= 1.003
%-------------------------
BK7glass= importdata("BK7glass.csv", ",",1)
nt=BK7glass.data(:,2)

% refractive index of material
nt=BK7glass.data(:,2)
%-------------------------
%angle of transmission(using snell)
tt= asin(ni.*cos(ti)./nt)
%-------------------------

rs=(ni.*cos(ti)-nt.*cos(tt))./(ni.*cos(ti)+nt.*cos(tt))
rp=(nt.*cos(tt)-ni.*cos(tt))./(nt.*cos(ti)+ni.*cos(tt))
ts=(2*ni.*cos(ti))./(ni.*cos(ti)+nt.*cos(ti))
tp=(2*ni.*cos(ti))./(nt.*cos(ti)+ni.*cos(ti))
%-------------------------
Rp=abs(rp).^2 %i tried to use a simple for loop to do this but matlab is complicated like that
Rs=abs(rs).^2
Tp=abs(tp).^2 
Ts=abs(ts).^2
%-------------------------
R=(Rp+Rs)/2
T=(Tp+Ts)/2
plt= plot(wl,R,wl,T)
xlim([0.29 1])
ylim([0 0.7])


