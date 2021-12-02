function Y = RungeKutta4( t,y0,F )

Lt=length(t);
h=t(2)-t(1);
Y=zeros(length(y0),Lt);
Y(:,1)=y0;
for j=1:Lt-1
   k1=h*F(t(j),Y(:,j));
   k2=h*F(t(j)+h/2,Y(:,j)+k1/2);
   k3=h*F(t(j)+h/2,Y(:,j)+k2/2);
   k4=h*F(t(j)+h,Y(:,j)+k3);
   Y(:,j+1)=Y(:,j)+(1/6)*(k1+2*k2+2*k3+k4);
end

end