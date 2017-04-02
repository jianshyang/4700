function [] = V_E()
clear;
close all;
gap=1/30;%the gap between potentials
W=1;
L=2;
%the ratio between the length and width

co=1;
ci=1e-2;
%set the conductivities

%normalize the conductivities so that inceasing mesh density will not
%chenge the total resistance over the rectangualar plate
c1=co.*gap;
c2=ci.*gap;

%boundary conditions
bc1=0.1;
bc2=0;

%locate conductivities over the plate
for j=1:1:W/gap+1
    for i=1:1:L/gap
        if i>0.4*L/gap && i<=0.6*L/gap+1
            if j<=0.4*W/gap+1 || j>0.6*W/gap+1
                C(j,i)=c2;
            else
                C(j,i)=c1;
            end
        else
             C(j,i)=c1;
        end
    end
end


[G,B]=FourFixedPoint(1,0,L/gap+1,W/gap+2,0,L/gap,1./C,1);
y2=G\B;
k=1;

for j=1:W/gap+1
   for i=1:L/gap+1
            V(j,i)=y2(k);
            k=k+1;
   end
end

[x1,y1]=meshgrid(0:L/gap,0:W/gap);
%calculate the E field
[ux,uy]=gradient(V,gap);

figure(1)
mesh(x1,y1,V);
title('Electrostatic Potential with a bottle-neck')
xlabel('x');
ylabel('y');
zlabel('Potential (V)');
hold off;

figure(2)
quiver(x1,y1,-ux,-uy,1.5);
xlim([0 L/gap]);
ylim([0 W/gap]);
title('Electric field  with a bottle-neck')
xlabel('x');
ylabel('y');
zlabel('E');
hold off;


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
                   G(h,h-1)=1/density(n,k-1);
                   G(h,h+1)=1/density(n,k);
                   G(h,h)=-(1/density(n,k-1)+1/density(n,k));
                   B(h,1)=0;
                elseif k<bnb 
                    B(h,1)=-bc1/density(n-1,k-1);
                    G(h,h-1)=1/density(n,k-1);
                    G(h,h+1)=1/density(n,k);
                    G(h,h)=-(1/density(n-1,k-1)+1/density(n,k-1)+1/density(n,k)+1/density(n+1,k-1));
                    G(h,h+(nx))=1/density(n+1,k-1); 
                elseif k>bne 
                    B(h,1)=-bc2/density(n-1,k-1);
                    G(h,h-1)=1/density(n,k-1);
                    G(h,h+1)=1/density(n,k);
                    G(h,h)=-(1/density(n-1,k-1)+1/density(n,k-1)+1/density(n,k)+1/density(n+1,k-1));
                    G(h,h+(nx))=1/density(n+1,k-1); 
                end
            elseif n==ny-1
                if k>=bnb && k<=bne
                   G(h,h-1)=1/density(n,k-1);
                   G(h,h+1)=1/density(n,k);
                   G(h,h)=-(1/density(n,k-1)+1/density(n,k));
                   B(h,1)=0;
                elseif k<bnb 
                    B(h,1)=-bc1/density(n+1,k-1);
                    G(h,h-1)=1/density(n,k-1);
                    G(h,h+1)=1/density(n,k);
                    G(h,h)=-(1/density(n-1,k-1)+1/density(n,k-1)+1/density(n,k)+1/density(n+1,k-1));
                    G(h,h-(nx))=1/density(n-1,k-1); 
                elseif k>bne 
                    B(h,1)=-bc2/density(n+1,k-1);
                    G(h,h-1)=1/density(n,k-1);
                    G(h,h+1)=1/density(n,k);
                    G(h,h)=-(1/density(n-1,k-1)+1/density(n,k-1)+1/density(n,k)+1/density(n+1,k-1));
                    G(h,h-(nx))=1/density(n-1,k-1);                 
                end
            else
                B(h,1)=0;
                G(h,h-1)=1/density(n,k-1);
                G(h,h+1)=1/density(n,k);
                G(h,h)=-(1/density(n,k-1)+1/density(n,k)+1/density(n-1,k-1)+1/density(n+1,k-1));
                G(h,h-nx)=1/density(n-1,k-1);
                G(h,h+(nx))=1/density(n+1,k-1); 
            end
            
    end
end
end
    