
function []=MFP_with_field()
close all;
clear;
clc;

k=1.381e-23;%boltzmann constant;
m=0.26*9.109e-31;%effective mass;
T=300;%temperature in K;
xbound=200e-9;
ybound=100e-9;
numofelectron=1000;
tau=0.2e-12;
timestep=2000;
vth=sqrt(2.*k.*T/m);%the thermal velocity in 2 degrees of freedom;
%mean free path
lamuda=vth*tau;
%set the fixed time interval that the electron can move less than 1/100 of the region size in each timestep
dt=ybound/200/vth;
clr=rand(20,3);
V=0.1;
q=1.60217662e-19;%charge
ec=10e15/(0.01*0.01);%electron concentration
Ex=V/xbound;%field in x-direction
Fx=-q.*Ex;%force in x-direction
ax=Fx/m;%acceleartion in x direction;
ay=0;

%uniformly set the random location of each electron
x=location(numofelectron,xbound);
y=location(numofelectron,ybound);

%use the maxwell-boltzmann distribution for each velocity conmponent with
%the average speed of vth
vx=direction(numofelectron,k,T,m);
vy=direction(numofelectron,k,T,m);

%current
I(1)=q.*ec.*abs(mean(vx)).*ybound;







figure(1);
hold on
title('Trajectories with Scattering');
xlabel('X (m)');
ylabel('Y (m)');
%set the region of the plot
xlim([0 200e-9]);
ylim([0 100e-9]);
%Array Tp is used to track the temperature of each timestep
%use display to show the temperature in the command for speeding up the
%simulation instead of using legend


for i=2:timestep
    %check if collisions occur. if Pscat > rand() then the particle
    %scatters. index save the indices of the electrons which need to be
    %re-thermalize
    vx(1,:)=vx(1,:)+ax.*dt;
    vy(1,:)=vy(1,:)+ay.*dt;
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
    
    %if a electron hit the y bound, the electron relfects by assign a equal
    %negative speed in y-direction. Determine how long the electron
    %excesses from the y-bound, then give the correct y-location 
    vy(find(tempy>100e-9|tempy<0)) =-vy(find(tempy>100e-9|tempy<0));
    tempy(find(tempy>100e-9))= (200e-9)-tempy(find(tempy>100e-9));
    tempy(find(tempy<0))= -tempy(find(tempy<0));
    
    %save the adjusted loction into the location arrays
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
    I(i)=q.*ec.*abs(mean(vx)).*ybound;
    pause(0.00001);

end



hold off;


figure(2)
plot((0:timestep-1).*dt,I);
title('Drift current vs. time');
xlabel('time (s)');
ylabel('Drift current (A)');

figure(3)
%plot the electron density map
histogram2(x(timestep,:),y(timestep,:),[50 50],'FaceColor','flat')
colorbar
xlim([0 200e-9]);
ylim([0 100e-9]);
title('The electron density map');
xlabel('X (m)');
ylabel('Y (m)');

figure(4)
[xq,yq]=meshgrid(0:0.1e-9:200e-9,0:0.1e-9:100e-9);
Z = griddata(x(timestep,:),y(timestep,:),(vx(1,:).^2+vy(1,:).^2).*m./k./2,xq,yq);
mesh(xq,yq,Z)
xlim([0 200e-9]);
ylim([0 100e-9]);
colorbar
title('The temperature map');
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
