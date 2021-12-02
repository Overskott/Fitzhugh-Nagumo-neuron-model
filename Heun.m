function Y = Heun( t,y0,F )

Lt=length(t);
h=t(2)-t(1);
Y=zeros(length(y0),Lt);
Y(:,1)=y0;
for j=1:Lt-1
   k1=F(t(j),Y(:,j));
   Yeuler=Y(:,j)+h*k1;
   k2=F(t(j+1),Yeuler);
   Y(:,j+1)=Y(:,j)+(h/2)*(k1+k2);
end

end

