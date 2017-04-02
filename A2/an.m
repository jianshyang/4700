function [] = an(  )
Vo=1;
temp=0;
a=20;
b=15;
[x,y]=meshgrid(-15:15,0:20);
f1=@(n) cosh((n.*pi.*abs(x))./a)./cosh((n.*pi.*b)./a);
f2=@(n) sin(n.*abs(y).*pi./a);
for n=1:2:51
    figure(1)
    temp=temp+4.*Vo./pi.*(1./n).*f1(n).*f2(n);
    mesh(x,y,temp);
    pause(1);
end

end

