%This script runs a first order Euler (FOE) integrator and a Verlet
%integrator using the same initial conditions and time step in order to
%show how they diverge. This script was created for a post about Verlet
%integration that will be posted on www.gereshes.com on 2018.07.06
%
%Ari Rubinsztejn
%ari@gereshes.com
%www.gereshes.com

close all
clear all
clc

%% Set these
dt=.1
planets=5;
finalTime=1000;% Time in seconds
G=.01%6.67408*(10^-11);%.01 Makes for cool plots, units exist but they arent important here


%% Initializing stuff
ICPos = rand(planets*2,1)*10;
ICVel = rand(planets*2,1)*.01;
dot=zeros(size(ICVel));
time=0:dt:finalTime;
timeSteps=length(time);
eulerPosHolder=zeros(2*planets,timeSteps);
eulerVelHolder=zeros(2*planets,timeSteps);
verletHolder=zeros(2*planets,timeSteps);
eulerPosHolder(:,1)=ICPos;
eulerVelHolder(:,1)=ICVel;
verletHolder(:,1)=ICPos;

%% Forward Euler step the equation
for c=2:timeSteps
    eulerPosHolder(:,c)=eulerPosHolder(:,c-1)+(eulerVelHolder(:,c-1)*dt);
    
    for d=1:2:2*planets
        acc=[0;0];
        pos=eulerPosHolder(d:d+1,c-1);
        for e=1:2:2*planets
            if(d==e)
            else
                posDiff=eulerPosHolder(e:e+1,c-1)-pos;
                acc=acc+(posDiff./(.001+((norm(posDiff))^3)));
            end
        end
        dot(d:d+1)=G*acc;
    end
    eulerVelHolder(:,c)=eulerVelHolder(:,c-1)+(dot*dt);
    c/timeSteps
end


%% Verlet integrate the equation forward
verletHolder(:,2)=eulerPosHolder(:,2);%you need two positions for 
for c=3:timeSteps
    for d=1:2:2*planets
        acc=[0;0];
        pos=verletHolder(d:d+1,c-1);
        for e=1:2:2*planets
            if(d==e)
            else
                posDiff=verletHolder(e:e+1,c-1)-pos;
                acc=acc+(posDiff./(.001+((norm(posDiff))^3)));
            end
        end
        dot(d:d+1)=G*acc;
    end
    verletHolder(:,c)=(2*verletHolder(:,c-1))-verletHolder(:,c-2)+(dot*dt*dt); %This is the verlet integrator step
    c/timeSteps +1
end

%% Plotting

%Plot the FOE
figure()
hold on
for c=1:planets
    plot(eulerPosHolder((2*c)-1,:),eulerPosHolder(2*c,:))
end
grid on
grid minor
title('First Order Euler (FOE) Integration')
xlabel('X Position')
ylabel('Y Position')

%Plot the verlet
figure()
hold on
for c=1:planets
    plot(verletHolder((2*c)-1,:),verletHolder(2*c,:))
end
grid on
grid minor
title('Verlet Integration')
xlabel('X Position')
ylabel('Y Position')

%Plot the difference 
figure()
hold on
for c=1:planets*2
    %plot(time,(verletHolder(c,:)-eulerPosHolder(c,:))*2./(verletHolder(c,:)+eulerPosHolder(c,:)))
    plot(time,(verletHolder(c,:)-eulerPosHolder(c,:)))
end
grid on
grid minor
title('Difference between Verlet and FOE Integration')
xlabel('Time (Seconds)')
ylabel('Difference')
