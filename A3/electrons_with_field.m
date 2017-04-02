function []=electrons_with_field()
close all;
clear;
clc;

k=1.381e-23;%boltzmann constant;
m=0.26*9.109e-31;%effective mass;
T=300;%temperature in K;
xbound=200e-9;
ybound=100e-9;
numofelectron=2000;
tau=0.2e-12;
timestep=2000;
vth=sqrt(2.*k.*T/m);%the thermal velocity in 2 degrees of freedom;
%mean free path
lamuda=vth*tau;
%set the fixed time interval that the electron can move less than 1/100 of the region size in each timestep
dt=ybound/200/vth;
clr=rand(20,3);
gap=0.1;

s=2;
if s==1
    bc=0.1;
elseif s==2
    bc=0.8;
end

q=1.60217662e-19;%charge
ec=10e15/(0.01*0.01);%electron concentration
[ux,uy]=V_E(gap,bc);%field
Ex=ux./1e-7;
Ey=uy./1e-7;
Fx=-q.*Ex;%force in x-direction
Fy=-q.*Ey;%force in y-direction
ax=Fx/m;%acceleartion in x direction;
ay=Fy/m;



%uniformly set the random location of each electron
x=location(numofelectron,xbound);
y=location(numofelectron,ybound);



%use a for loop, a while loop and a if statement to check if any electron
%is in one of boxes region
for j=1:numofelectron
    check=true;
    while check
        %while loop keeps checking if a electronis in one of boxes region.
        %if it is, re-assign location and check again. if no, start
        %checking the next one
        if (x(1,j)>80e-9 && x(1,j)<120e-9) && ((y(1,j)>60e-9)||(y(1,j)<40e-9))
            x(1,j)=location(1,xbound);
            y(1,j)=location(1,ybound);
        else
            check = false;
        end
    end
end



%use the maxwell-boltzmann distribution for each velocity conmponent with
%the average speed of vth
vx=direction(numofelectron,k,T,m);
vy=direction(numofelectron,k,T,m);





figure(1);
hold on
title('Trajectories with Scattering and two boxse');
xlabel('X (m)');
ylabel('Y (m)');
%locate two boxes on the plot
rectangle('Position',[80e-9,60e-9,40e-9,40e-9]);
rectangle('Position',[80e-9,0,40e-9,40e-9]);

%Array Tp is used to track the temperature of each timestep
%use display to show the temperature in the command for speeding up the
%simulation instead of using legend
xlim([0 200e-9]);
ylim([0 100e-9]);





for i=2:timestep
    %check if collisions occur. if Pscat > rand() then the particle
    %scatters. index save the indices of the electrons which need to be
    %re-thermalize
    vx(1,:)=vx(1,:)+ax(sub2ind(size(ax),round((y(i-1,:)-rem(y(i-1,:),gap.*1e-7)).*1e7./gap+1),round((x(i-1,:)-rem(x(i-1,:),gap.*1e-7)).*1e7./gap+1))).*dt;
    vy(1,:)=vy(1,:)+ay(sub2ind(size(ax),round((y(i-1,:)-rem(y(i-1,:),gap.*1e-7)).*1e7./gap+1),round((x(i-1,:)-rem(x(i-1,:),gap.*1e-7)).*1e7./gap+1))).*dt;
    
    index=find((1-exp(-dt./tau))>rand(1,numofelectron));
    [a b]=size(index);
    %re-thermalize
    vx(1,index)=direction(b,k,T,m);
    vy(1,index)=direction(b,k,T,m); 
    %use two temp arrays to save the new location 
    tempx=x(i-1,:)+vx(1,:).*dt;
    tempy=y(i-1,:)+vy(1,:).*dt;
    %if a electron hit the x bound, the electron jumps to the opposite edge
    tempx(find(tempx>200e-9))= tempx(find(tempx>200e-9))- 200e-9;
    tempx(find(tempx<0))= (200e-9)+tempx(find(tempx<0));
    
    %this for loop checks if a electron hits a box. If hits, reflect the
    %electron
    for j=1:numofelectron
        
        if tempx(j)>=80e-9 && tempx(j)<=85e-9
            if tempy(j)>=60e-9 || tempy(j) <= 40e-9
                vx(1,j)=-vx(1,j);
                tempx(j)=160e-9-tempx(j);
            end
        elseif tempx(j)>=115e-9 && tempx(j)<=120e-9
            if tempy(j)>=60e-9 || tempy(j) <= 40e-9
                vx(1,j)=-vx(1,j);
                tempx(j)=240e-9-tempx(j);
            end    
        end
        
        if tempy(j)>=60e-9 && tempy(j) <=65e-9
            if tempx(j)>=80e-9 && tempx(j) <= 120e-9
                vy(1,j)=-vy(1,j);
                tempy(j)=120e-9-tempy(j);
            end
        elseif tempy(j)>=35e-9 && tempy(j) <=40e-9
            if tempx(j)>=80e-9 && tempx(j) <= 120e-9
                vy(1,j)=-vy(1,j);
                tempy(j)=80e-9-tempy(j);
            end
        end
                   
            
   
    end


    %if a electron hit the y bound, the electron relfects by assign a equal
    %negative speed in y-direction. Determine how long the electron
    %excesses from the y-bound, then give the correct y-location 
    vy(find(tempy>100e-9|tempy<0)) =-vy(find(tempy>100e-9|tempy<0));
    tempy(find(tempy>100e-9))= (200e-9)-tempy(find(tempy>100e-9));
    tempy(find(tempy<0))= -tempy(find(tempy<0));
    x(i,:)=tempx;
    y(i,:)=tempy;

    %plot first 5 electrons
    for  p=1:5
        %Because of the x boundary condition, plot function may a line
        %across the the region. This if statement prevents that happens
        if abs(x(i-1,p)-x(i,p)) < 100e-9
            plot(x(i-1:i,p),y(i-1:i,p),'Color',clr(p,:),'LineWidth',1);
        end     
    end

    pause(0.0001);

end
    
hold off;


figure(2)
%plot the electron density map
histogram2(x(timestep,:),y(timestep,:),[50 50],'FaceColor','flat')
colorbar
xlim([0 200e-9]);
ylim([0 100e-9]);
title('The electron density map');
xlabel('X (m)');
ylabel('Y (m)');








end

function [loc]= location(Numofpoint,bound)
%uniformly generate random x or y loction
loc=bound.*(rand(1,Numofpoint));
end
function [dir]=direction(Numofpoint,k,T,m)
%use the maxwell-boltzmann distribution for each velocity conmponent
dir=randn(1,Numofpoint).*sqrt(k.*T./m);
end

function [Ex,Ey] = V_E(gap,bc)


W=1;
L=2;
co=1;
ci=1e-2;
%normalize the conductivities so that inceasing mesh density will not
%chenge the total resistance over the rectangualar plate
c1=co.*gap;
c2=ci.*gap;

bc1=bc;
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


[G,B]=FourFixedPoint(bc1,bc2,L/gap+1,W/gap+2,0,L/gap,1./C,1);
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

Ex=-ux;
Ey=-uy;



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
