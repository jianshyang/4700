function [  ] = laplace(  )
L=300;
W=200;
size=1;
c1=1;
c2=10e-2;
for j=1:size:W+1
    for i=1:size:L+1
        if i>=L./2-L./10 && i<L./2+L/10
            if j<=W./2-W./10 || j> W./2+W./10 
                C(j,i)=c2;
            else
                C(j,i)=c1;
            end
        else
             C(j,i)=c1;
        end
    end
end
[x,y]=meshgrid(0:L,0:W);
figure(1);
mesh(x,y,C)

end
    