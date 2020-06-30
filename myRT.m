function[R,T,Rp,Rs,Tp,Ts]=myRT(ti,Ni,Nt)
%-----------fresnel--------------
rp=(((Nt/Ni).^2)*cos(ti)-((Nt/Ni).^2-sin(ti)^2).^0.5)./(((Nt/Ni).^2)*cos(ti)+((Nt/Ni).^2+sin(ti)^2).^0.5)
rs=(cos(ti)-((Nt/Ni).^2-sin(ti)^2).^0.5)./(cos(ti)+((Nt/Ni).^2+sin(ti)^2).^0.5)
tp=2*cos(ti)./(cos(ti)+((Nt/Ni).^2-sin(ti)^2).^0.5)
ts=2*cos(ti)./((Nt/Ni).*cos(ti)+(1-((Nt/Ni).*sin(ti)).^2).^0.5)
%-----------polarization----------
Rp=abs(rp).^2 %i tried to use a simple for loop to do this but matlab is complicated like that
Rs=abs(rs).^2
Tp=(((Nt/Ni).^2+sin(ti)^2).^0.5).*(abs(tp).^2)./cos(ti) 
Ts=(((Nt/Ni).^2+sin(ti)^2).^0.5).*(abs(ts).^2)./cos(ti) 
%------------unpolarized----------
R=(Rp+Rs)/2
T=(Tp+Ts)/2
end