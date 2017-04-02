function [ ] = laplace(  )
clear;
clc;
close all;

%solve a 1-D case when L=30 and conductivity = 1 per length 
%the boundary condition is V=1 at x=0,and V=0 at x=L
[G,B]=TwoFixedPoint(1,0,31,1,1);
y=G\B;
x=0:1:30;
figure(1)
plot(x,y)
title('Electrostatic Potential vs. Length')
xlabel('Length');
ylabel('Potential (V)');


%solve a 2-D case when L=30, W=20 and conductivity = 1 per length 
%the boundary condition is V=1 at x=0 and x=L, and V=0 at y=0 and y=W
[A B]=FourFixedPoint(1,1,0,0,31,20,1,1);
y2=A\B;
k=1;
%sort the potentials into a matrix for generating a mesh plot
for j=1:21
   for i=1:31
       if j ==1 || j==21
            V(j,i)=0;
            if i==1 || i==31
                V(j,i)=1;
            end
       else
            V(j,i)=y2(k);
            k=k+1;
       end
   end
end

[x1,y1]=meshgrid(0:30,0:20);



figure(2)
mesh(x1,y1,V);
title('Electrostatic Potential over a 3/2 Rectangular region')
xlabel('x');
ylabel('y');
zlabel('Potential (V)')


%solve 2D Case analytically
Vo=1;%boundard condition
temp=0;
a=20;%W=20
b=15;%L=30 from -15 to 15
[x,y]=meshgrid(-15:15,0:20);
%creat functions of each term in the analytical equation for simplicity
f1=@(n) cosh((n.*pi.*abs(x))./a)./cosh((n.*pi.*b)./a);
f2=@(n) sin(n.*abs(y).*pi./a);
figure(3)
%start the iteration and figure out it should stop
for n=1:2:301
    temp=temp+4.*Vo./pi.*(1./n).*f1(n).*f2(n);
    mesh(x,y,temp);
    xlim([-15 15]);
    ylim([0 20]);
    title('analytical solution of Electrostatic Potential over a 3/2 Rectangular region')
    xlabel('x');
    ylabel('y');
    zlabel('Potential (V)')
    %track the error between the numerical solution and analytical solution
    error=V-temp;
    pause(0.001);
end
%return the final error
error=V-temp;

end


function [G,B]=TwoFixedPoint(bc1,bc2,n,density,stepsize)
%{
bc1,bc2: boundary condition
n:number of potiential points 
densisty:resistance densisty,ohm/step
stepsize: the distance between two closet points
%}
for k=1:stepsize:n
    if k==1
        %when x=0, set its boundary condition
        B(k,1)=bc1;
        G(k,1)=1;
    elseif k==n
        %when x=L, set its boundary condition
        B(k,1)=bc2;
        G(k,k)=1;
    else
        B(k,1)=0;
        G(k,k-1)=1/density;
        G(k,k)=-2/density;
        G(k,k+1)=1/density;  
    end
end

end

function [G,B]=FourFixedPoint(bc1,bc2,bc3,bc4,nx,ny,density,size)
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
                %when x=0, set its boundary condition
                G(h,h)=1;
                B(h,1)=bc1;
            elseif k==nx
                %when x=L, set its boundary condition
                G(h,h)=1;
                B(h,1)=bc2;
            elseif n==1
                %set the equation when the point is 1 unit above the lower
                %y-bound
                B(h,1)=-bc4/density;
                G(h,h-1)=1/density;
                G(h,h+1)=1/density;
                G(h,h)=-4/density;
                G(h,h+(nx))=1/density; 
            elseif n==ny-1
                %set the equation when the point is 1 unit below the upper
                %y-bound
                B(h,1)=-bc3/density;
                G(h,h+1)=1/density;
                G(h,h-1)=1/density;
                G(h,h)=-4/density;
                G(h,h-(nx))=1/density; 
            else
                B(h,1)=0;
                G(h,h+1)=1/density;
                G(h,h-1)=1/density;
                G(h,h)=-4/density;
                G(h,h-nx)=1/density;
                G(h,h+nx)=1/density; 
            end
            
    end
end
        
            
            
        
   
        
        
    

end

