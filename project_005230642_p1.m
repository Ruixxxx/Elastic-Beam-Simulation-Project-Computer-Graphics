%%Rigid sphere falling in viscous flow
%Use the implicit and explicit Euler method to calculate the position and
%velocity of the sphere in viscous fluid vs. time with initial conditions
%and plot them. Analyse the accuary of these two simulations by computing
%the theoretical terminal velocity and the position of the sphere when the
%viscosity is 0 and comparing them with the results from simulations.

%%Typically scripts begin with
%{
Rui Xu
005230642
December 14,2018
%}
clc;clear all;close all;

%%Assign the variables
%Attribute variables
mu=1000;        %viscosity of the fluid
rho_fluid=1000; %density of the fluid
rho_metal=7000; %density of the sphere
R=0.035;        %radius of the sphere
[W,m]=WeightMass(R,rho_metal,rho_fluid); %weight and mass

% tic
%%Implicit simulation
%Time counter
dt=0.01;
tfinal=1;
%Calculate the steps in for-loop
tstep=ceil(tfinal/dt);

%Initialize the position and velocity of the sphere
yc=zeros(1,tstep+1);
yc_dot=zeros(1,tstep+1);

%Initial conditions
yc(1)=0;
yc_dot(1)=0;

%Calculate the position and velocity of the sphere implicitly
for k=1:tstep
    %Initialize the error
    error=10;
    %Initialize the tolerance
    tol=1e-3;
    %Set the starting point
    q=yc(k);
    while error>tol
        %Compute the function
        f=m/dt*((q-yc(k))/dt-yc_dot(k))+W+(6*pi*mu*R)*(q-yc(k))/dt;
        %Compute the Jacobian
        J=m/(dt^2)+(6*pi*mu*R)/dt;
        %Newton's update
        dq=J\f;
        q=q-dq;
        %Compute the error
        error=abs(f);
    end
        %Position
        yc(k+1)=q;
        %Velocity
        yc_dot(k+1)=(yc(k+1)-yc(k))/dt;
end

%Plot the position and velocity vs. time for implicit simulation
t=0:dt:tfinal;
%Create the plot
figure(1)
plot(t,yc);
%Format the plot
xlabel('time (sec)');
ylabel('position (m)');
title('the Position vs. Time for Implicit Simulation');
%Create the plot
figure(2)
plot(t,yc_dot);
%Format the plot
xlabel('time (sec)');
ylabel('velocity (m/sec)');
title('the Velocity vs. Time for Implicit Simulation');
% toc

% tic
%%Explicit simulation
%Time counter
dT=0.0001;
Tfinal=1;
%Calculate the steps in for-loop
Tstep=ceil(Tfinal/dT);

%Initialize the position and velocity of the sphere
YC=zeros(1,Tstep+1);
YC_dot=zeros(1,Tstep+1);

%%Initial conditions
YC(1)=0;
YC_dot(1)=0;

%%Calculate the position and velocity if the sphere explicitly
for k=1:Tstep
        %Position
        YC(k+1)=((dT/m)*(-W-6*pi*mu*R*YC_dot(k))+YC_dot(k))*dT+YC(k);
        %Velocity
        YC_dot(k+1)=(YC(k+1)-YC(k))/dT;
end

%%Plot the position and velocity vs. time for explicit simulation
%Create the plot
T=0:dT:Tfinal;
figure(3)
plot(T,YC);
%Format the plot
xlabel('time (sec)');
ylabel('position (m)');
title('the Position vs. Time for Explicit Simulation');
%Create the plot
figure(4)
plot(T,YC_dot);
%Format the plot
xlabel('time (sec)');
ylabel('velocity (m/sec)');
title('the Velocity vs. Time for Explicit Simulation');
% toc

% %%Calculte the theoretical terminal velocity
% velocity_terminal=-W/(6*pi*mu*R);
% %Display the terminal velocities
% fprintf('The terminal velocity from implicit simulation is %f.\n',yc_dot(tstep+1));
% fprintf('The terminal velocity from explicit simulation is %f.\n',YC_dot(Tstep+1));
% fprintf('The theoretical terminal velocity is %f.\n',velocity_terminal);
% %Calculate the position of the freefall sphere
% tt=0:0.01:1;
% g=9.8;
% yc_freefall=-1/2*g*tt.^2;
% %%Plot the position of the free fall sphere
% %Create the plot
% figure(5)
% plot(tt,yc_freefall);
% %Format the plot
% xlabel('time (sec)');
% ylabel('position (m)');
% title('the Position vs. Time for Freefall');