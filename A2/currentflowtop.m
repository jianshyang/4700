function[] = currentflowtop()
clear;
clc;
close all;
%calculate the current
I1=currentflow(30,20,15,15,1,1e-2)
%manually select each investigation case
i=0;
if i==1
    %mesh density
    close all;
    I2=currentflow(60,40,30,30,1,1e-2)
    close all;
    I3=currentflow(90,60,45,45,1,1e-2)
    close all;
    figure(5)
    plot([30*20 60*40 90*60],[I1 I2 I3])
    title('Current vs Mesh Density')
    xlabel('mesh density');
    ylabel('Current');

elseif i==2
    %size of bottle-neck
    close all;
    I4=currentflow(30,20,14,16,1,1e-2)
    close all;
    I5=currentflow(30,20,13,17,1,1e-2)
    close all;
    I6=currentflow(30,20,12,18,1,1e-2)
    close all;
    I7=currentflow(30,20,11,19,1,1e-2)
    close all;
    I8=currentflow(30,20,10,20,1,1e-2)
    close all;
    I9=currentflow(30,20,9,21,1,1e-2)
    close all;
    figure(6)
    plot([2 4 6 8 10 12 14],[I1 I4 I5 I6 I7 I8 I9])
    title('Current vs bottle-neck size')
    xlabel('bottle-neck size');
    ylabel('Current');
    xlim([2 14])
elseif i==3
    %varing the conductivity
    close all;
    I13=currentflow(30,20,15,15,1,10e-2)
    close all;
    I10=currentflow(30,20,15,15,1,5e-2)
    close all;
    I11=currentflow(30,20,15,15,1,5e-3)
    close all;
    I12=currentflow(30,20,15,15,1,1e-3)
    close all;
    figure(7)
    plot([10e-2 5e-2 1e-2 5e-3 1e-3],[I13 I10 I1 I11 I12])
    title('Current vs conductivity')
    xlabel('conductivity');
    ylabel('Current');
end

end
function [current] = currentflow(L,W,bnb,bne,co,ci)

size=1;

%normalize the conductivities so that inceasing mesh density will not
%chenge the total resistance over the rectangualar plate
c1=co./(L./30);
c2=ci./(L./30);

bc1=1;
bc2=0;
%locate conductivities over the plate
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
title('Conductivity ver a 3/2 Rectangular region')
xlabel('x');
ylabel('y');
zlabel('Conductivity (ohm^-1)')
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
%assume the total current flow is the amount of current flow in/out of the
%bottle-neck region
I=sum((V(2:W,bnb-1)-V(2:W,bnb)).*C(2:W,bnb-1));
[x1,y1]=meshgrid(0:L,0:W);

%calculate the E field
[ux,uy]=gradient(V(2:W,:),0.1);

figure(2)
mesh(x1,y1,V);
title('Electrostatic Potential over a 3/2 Rectangular region with a bottle-neck')
xlabel('x');
ylabel('y');
zlabel('Potential (V)');

figure(3)
hold on;
title('Electric field over a 3/2 Rectangular region with a bottle-neck')
xlabel('x');
ylabel('y');
zlabel('E');
quiver(x1(2:W,:),y1(2:W,:),-ux,-uy,1.3);
rectangle('Position',[0,0,L,W]);
xlim([0 L]);
ylim([0 W]);
hold off;

figure(4)
hold on
title('Current density over a 3/2 Rectangular region with a bottle-neck')
xlabel('x');
ylabel('y');
zlabel('J');
%calculate the current density 
quiver(x1(2:W,:),y1(2:W,:),-ux.*C(2:W,:),-uy.*C(2:W,:),1.3);
rectangle('Position',[0,0,L,W]);
xlim([0 L]);
ylim([0 W]);
hold off;

current = I;
end

function [G,B]=FourFixedPoint(bc1,bc2,nx,ny,bnb,bne,density,size)
%{
bc1:boundary condition at the left bound
bc2:boundary condition at the right bound
bc3:boundary condition at the top bound
bc4:boundary condition at the bottom bound
nx:number of points along x-axis,inclusive two bounds 
ny;number of points along y-axis,exclusive two bounds
densisty:resistance densisty,ohm/length
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
    