function[] = currentflowtop()
end
function [current ] = currentflow(L,W,bnb,bne,co,ci)
clear;
clc;
close all;
size=1;
c1=co./(L./30);
c2=ci./(L./30);
bc1=1;
bc2=0;
for j=1:size:W+1
    for i=1:size:L+1
        if i>L./2-L./10 && i<=L./2+L/10
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
mesh(x,y,C);
[G,B]=FourFixedPoint(1,0,L+1,W,bnb,bne,1./C,1);
y2=G\B;
k=1;

for j=1:W+1
   for i=1:L+1
       if j ==1 || j==W+1
            V(j,i)=0;
            if i<bnb
                V(j,i)=1;
            end
       else
            V(j,i)=y2(k);
            k=k+1;
       end
   end
end
I=sum((V(2:W,bnb-1)-V(2:W,bnb)).*C(2:W,bnb-1))
[x1,y1]=meshgrid(0:L,0:W);
[ux,uy]=gradient(V(2:W,:),0.1);
figure(2)
mesh(x1,y1,V);

figure(3)
hold on;

quiver(x1(2:W,:),y1(2:W,:),-ux,-uy,1.3);
rectangle('Position',[0,0,L,W]);
xlim([0 L]);
ylim([0 W]);
hold off;

figure(4)
hold on
quiver(x1(2:W,:),y1(2:W,:),-ux.*C(2:W,:),-uy.*C(2:W,:),1.3);
rectangle('Position',[0,0,L,W]);
xlim([0 L]);
ylim([0 W]);
hold off;

current = I;
end

function [G,B]=FourFixedPoint(bc1,bc2,nx,ny,bnb,bne,density,size)
%{
function [G,B]=FourFixedPoint(bc1,bc2,bc3,bc4,L,W,density,size)
bc1:boundary condition at the left bound
bc2:boundary condition at the right bound
bc3:boundary condition at the top bound
bc4:boundary condition at the bottom bound
L:the length
densisty:resistance densisty,ohm/step
stepsize: the distance between two closet points
%}
for n=1:size:ny-1
    for k=1:size:nx
            h=k+(n-1).*(nx);
            if k==1 
                G(h,h)=1;
                B(h,1)=bc1;
            elseif k==nx
                G(h,h)=1;
                B(h,1)=bc2;
            elseif n==1
                if k>=bnb && k<=bne
                   G(h,h-1)=1/density(n+1,k-1);
                   G(h,h+1)=1/density(n+1,k);
                   G(h,h)=-(1/density(n+1,k-1)+1/density(n+1,k));
                   B(h,1)=0;
                elseif k<bnb 
                    B(h,1)=-bc1/density(n,k);
                    G(h,h-1)=1/density(n+1,k-1);
                    G(h,h+1)=1/density(n+1,k);
                    G(h,h)=-(1/density(n,k)+1/density(n+1,k-1)+1/density(n+1,k)+1/density(n+2,k));
                    G(h,h+(nx))=1/density(n+2,k); 
                elseif k>bne 
                    B(h,1)=-bc2/density(n,k);
                    G(h,h-1)=1/density(n+1,k-1);
                    G(h,h+1)=1/density(n+1,k);
                    G(h,h)=-(1/density(n,k)+1/density(n+1,k-1)+1/density(n+1,k)+1/density(n+2,k));
                    G(h,h+(nx))=1/density(n+2,k); 
                end
            elseif n==ny-1
                if k>=bnb && k<=bne
                   G(h,h-1)=1/density(n+1,k-1);
                   G(h,h+1)=1/density(n+1,k);
                   G(h,h)=-(1/density(n+1,k-1)+1/density(n+1,k));
                   B(h,1)=0;
                elseif k<bnb 
                    B(h,1)=-bc1/density(n+2,k);
                    G(h,h-1)=1/density(n+1,k-1);
                    G(h,h+1)=1/density(n+1,k);
                    G(h,h)=-(1/density(n,k)+1/density(n+1,k-1)+1/density(n+1,k)+1/density(n+2,k));
                    G(h,h-(nx))=1/density(n,k); 
                elseif k>bne 
                    B(h,1)=-bc2/density(n+2,k);
                    G(h,h-1)=1/density(n+1,k-1);
                    G(h,h+1)=1/density(n+1,k);
                    G(h,h)=-(1/density(n,k)+1/density(n+1,k-1)+1/density(n+1,k)+1/density(n+2,k));
                    G(h,h-(nx))=1/density(n,k);                 
                end
            else
                B(h,1)=0;
                G(h,h-1)=1/density(n+1,k-1);
                G(h,h+1)=1/density(n+1,k);
                G(h,h)=-(1/density(n,k)+1/density(n+1,k-1)+1/density(n+1,k)+1/density(n+2,k));
                G(h,h-nx)=1/density(n,k);
                G(h,h+(nx))=1/density(n+2,k); 
            end
            
    end
end
end
    