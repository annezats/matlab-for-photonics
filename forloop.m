

rp=0:5
rs=0
tp=0
ts=0

fresnel={rp,rs,tp,ts}%this works only if vars are already assigned
big={'Rp','Rs','Tp','Ts'} %if vars to be assigned later, make them strings
disp(Rp)
for i=1:4
    big{i}=abs(fresnel{i}).^2 %{} access info inside cell, () access cell itself
end
disp(Ts)