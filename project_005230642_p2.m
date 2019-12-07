%%Rigid spheres and elastic beam falling in viscous flow
%Use the implicit and explicit Euler method to calculate the positions and
%velocities of the spheres in viscous fluid vs. time with initial
%conditions. Plot the position and velocity of the middle node vs. time
%and the turning angle vs. time.

%%Typically scripts begin with
%{
Rui Xu
005230642
December 14,2018
%}
clc;clear all;close all;

%%Assign the variables
r0=0.001;       %cross-sectional radius of the beam
E=1e9;          %Young's modulus
EA=E*pi*r0^2;   %stretching stiffness
EI=E*pi*r0^4/4; %bending stiffness
l=0.1;          %length of the beam
l_k=l/2;        %length of each segment between two nodes

mu=1000;        %viscosity of the fluid
rho_fluid=1000; %density of the fluid
rho_metal=7000; %density of the sphere
R1=0.005;       %radii of the spheres
R2=0.025;
R3=0.005;

[W1,m1]=WeightMass(R1,rho_metal,rho_fluid); %weights and masses
[W2,m2]=WeightMass(R2,rho_metal,rho_fluid);
[W3,m3]=WeightMass(R3,rho_metal,rho_fluid);

% tic
%%Implicit simulation
%Time counter
dt=0.01;
tfinal=10;
%Calculate the steps in for-loop
tstep=ceil(tfinal/dt);

%Initialize the positions and velocities of the spheres
x1=zeros(1,tstep+1);x1_dot=zeros(1,tstep+1);
y1=zeros(1,tstep+1);y1_dot=zeros(1,tstep+1);
x2=zeros(1,tstep+1);x2_dot=zeros(1,tstep+1);
y2=zeros(1,tstep+1);y2_dot=zeros(1,tstep+1);
x3=zeros(1,tstep+1);x3_dot=zeros(1,tstep+1);
y3=zeros(1,tstep+1);y3_dot=zeros(1,tstep+1);
%Initial conditions
x1(1)=0;x1_dot(1)=0;
y1(1)=0;y1_dot(1)=0;
x2(1)=l_k;x2_dot(1)=0;
y2(1)=0;y2_dot(1)=0;
x3(1)=2*l_k;x3_dot(1)=0;
y3(1)=0;y3_dot(1)=0;

%Calculate the positions and velocities of the spheres implicitly
for k=1:tstep
    %Initialize the error
    error=10;
    %Initialize the tolerance
    tol=1e-3;
    %Set the starting points
    q=[x1(k);y1(k);x2(k);y2(k);x3(k);y3(k)];
    while error>tol
        %Initialize the functions
        f=zeros(6,1);
        Fs1=zeros(4,1);
        Fs2=zeros(4,1);
        Fb=zeros(6,1);
        %Initialize the Jacobian
        J=zeros(6,6);
        Js1=zeros(4,4);
        Js2=zeros(4,4);
        Jb=zeros(6,6);
        %Compute the functions
        Fs1=gradEs(q(1),q(2),q(3),q(4),l_k,EA);
        Fs2=gradEs(q(3),q(4),q(5),q(6),l_k,EA);
        Fb=gradEb(q(1),q(2),q(3),q(4),q(5),q(6),l_k,EI);
        f(1)=m1/dt*((q(1)-x1(k))/dt-x1_dot(k))+...
            (6*pi*mu*R1)*((q(1)-x1(k))/dt);
        f(2)=m1/dt*((q(2)-y1(k))/dt-y1_dot(k))+...
            W1+(6*pi*mu*R1)*((q(2)-y1(k))/dt);
        f(3)=m2/dt*((q(3)-x2(k))/dt-x2_dot(k))+...
            (6*pi*mu*R2)*((q(3)-x2(k))/dt);
        f(4)=m2/dt*((q(4)-y2(k))/dt-y2_dot(k))+...
            W2+(6*pi*mu*R2)*((q(4)-y2(k))/dt);
        f(5)=m3/dt*((q(5)-x3(k))/dt-x3_dot(k))+...
            (6*pi*mu*R3)*((q(5)-x3(k))/dt);
        f(6)=m3/dt*((q(6)-y3(k))/dt-y3_dot(k))+...
            W3+(6*pi*mu*R3)*((q(6)-y3(k))/dt);
        f(1:4,1)=f(1:4,1)+Fs1;
        f(3:6,1)=f(3:6,1)+Fs2;
        f=f+Fb;
        %Compute the Jacobian
        Js1=hessEs(q(1),q(2),q(3),q(4),l_k,EA);
        Js2=hessEs(q(3),q(4),q(5),q(6),l_k,EA);
        Jb=hessEb(q(1),q(2),q(3),q(4),q(5),q(6),l_k,EI);
        J(1,1)=m1/(dt^2)+(6*pi*mu*R1)/dt;
        J(2,2)=m1/(dt^2)+(6*pi*mu*R1)/dt;
        J(3,3)=m2/(dt^2)+(6*pi*mu*R2)/dt;
        J(4,4)=m2/(dt^2)+(6*pi*mu*R2)/dt;
        J(5,5)=m3/(dt^2)+(6*pi*mu*R3)/dt;
        J(6,6)=m3/(dt^2)+(6*pi*mu*R3)/dt;
        J(1:4,1:4)=J(1:4,1:4)+Js1';
        J(3:6,3:6)=J(3:6,3:6)+Js2';
        J=J+Jb';
        %Newton's update
        dq=J\f;
        q=q-dq;
        %Compute the error
        error=sum(abs(f));
    end
    %Position
    x1(k+1)=q(1);
    y1(k+1)=q(2);
    x2(k+1)=q(3);
    y2(k+1)=q(4);
    x3(k+1)=q(5);
    y3(k+1)=q(6);
    %Velocity
    x1_dot(k+1)=(x1(k+1)-x1(k))/dt;
    y1_dot(k+1)=(y1(k+1)-y1(k))/dt;
    x2_dot(k+1)=(x2(k+1)-x2(k))/dt;
    y2_dot(k+1)=(y2(k+1)-y2(k))/dt;
    x3_dot(k+1)=(x3(k+1)-x3(k))/dt;
    y3_dot(k+1)=(y3(k+1)-y3(k))/dt;
end

%%Plot the position and velocity of the middle node vs. time
t=0:dt:tfinal;
%Create the plot
figure(1)
plot(t,x2);
hold on
plot(t,y2);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('position (m)');
legend('x-position','y-position','Location','Best');
title('the Middle Node Position for Implicit Simulation');
%Create the plot
figure(2)
plot(t,x2_dot);
hold on
plot(t,y2_dot);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('velocity (m/sec)');
ylim([-9e-3,1e-3]);
legend('x-velocity','y-velocity','Location','Best');
title('the Middle Node Velocity for Implicit Simulation');

%%Plot the turning angle
%Calculate the turning angle
theta=acos((2*l_k^2-((2*x2-x1-x3).^2+(2*y2-y1-y3).^2))/(2*l_k^2));
%Create the plot
figure(3)
plot(t,theta);
%Format the plot
xlabel('time (sec)');
ylabel('angle (rad)');
ylim([0,1.2]);
title('the Turning Angle vs. Time for Implicit Simulation');
% toc

% tic
%%Explicit simulation
%Time counter
dT=1e-5;
Tfinal=10;
%Calculate the steps in for-loop
Tstep=ceil(Tfinal/dT);

%Initialize the positions and velocities of the sphere
X1=zeros(1,Tstep+1);X1_dot=zeros(1,Tstep+1);
Y1=zeros(1,Tstep+1);Y1_dot=zeros(1,Tstep+1);
X2=zeros(1,Tstep+1);X2_dot=zeros(1,Tstep+1);
Y2=zeros(1,Tstep+1);Y2_dot=zeros(1,Tstep+1);
X3=zeros(1,Tstep+1);X3_dot=zeros(1,Tstep+1);
Y3=zeros(1,Tstep+1);Y3_dot=zeros(1,Tstep+1);
%Initial conditions
X1(1)=0;X1_dot(1)=0;
Y1(1)=0;Y1_dot(1)=0;
X2(1)=l_k;X2_dot(1)=0;
Y2(1)=0;Y2_dot(1)=0;
X3(1)=2*l_k;X3_dot(1)=0;
Y3(1)=0;Y3_dot(1)=0;

%Calculate the positions and velocities of the spheres explicitly
for k=1:Tstep
    Fs1=gradEs(X1(k),Y1(k),X2(k),Y2(k),l_k,EA);
    Fs2=gradEs(X2(k),Y2(k),X3(k),Y3(k),l_k,EA);
    Fb=gradEb(X1(k),Y1(k),X2(k),Y2(k),X3(k),Y3(k),l_k,EI);
    %Position
    X1(k+1)=((dT/m1)*(-6*pi*mu*R1*X1_dot(k)-Fs1(1)-Fb(1))+X1_dot(k))*dT+X1(k);
    Y1(k+1)=((dT/m1)*(-W1-6*pi*mu*R1*Y1_dot(k)-Fs1(2)-Fb(2))+Y1_dot(k))*dT+Y1(k);
    X2(k+1)=((dT/m2)*(-6*pi*mu*R2*X2_dot(k)-Fs1(3)-Fs2(1)-Fb(3))+X2_dot(k))*dT+X2(k);
    Y2(k+1)=((dT/m2)*(-W2-6*pi*mu*R2*Y2_dot(k)-Fs1(4)-Fs2(2)-Fb(4))+Y2_dot(k))*dT+Y2(k);
    X3(k+1)=((dT/m3)*(-6*pi*mu*R3*X3_dot(k)-Fs2(3)-Fb(5))+X3_dot(k))*dT+X3(k);
    Y3(k+1)=((dT/m3)*(-W3-6*pi*mu*R3*Y3_dot(k)-Fs2(4)-Fb(6))+Y3_dot(k))*dT+Y3(k);
    %Velocity
    X1_dot(k+1)=(X1(k+1)-X1(k))/dT;
    Y1_dot(k+1)=(Y1(k+1)-Y1(k))/dT;
    X2_dot(k+1)=(X2(k+1)-X2(k))/dT;
    Y2_dot(k+1)=(Y2(k+1)-Y2(k))/dT;
    X3_dot(k+1)=(X3(k+1)-X3(k))/dT;
    Y3_dot(k+1)=(Y3(k+1)-Y3(k))/dT;
end

%%Plot the position and velocity of the middle node vs. time
T=0:dT:Tfinal;
%Create the plot
figure(4)
plot(T,X2);
hold on
plot(T,Y2);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('position (m)');
legend('x-position','y-position','Location','Best');
title('the Middle Node Position for Explicit Simulation');
%Create the plot
figure(5)
plot(T,X2_dot);
hold on
plot(T,Y2_dot);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('velocity (m/sec)');
ylim([-9e-3,1e-3]);
legend('x-velocity','y-velocity','Location','Best');
title('the Middle Node Velocity for Explicit Simulation');

%%Plot the turning angle
%Calculate the turning angle
THETA=acos((2*l_k^2-((2*X2-X1-X3).^2+(2*Y2-Y1-Y3).^2))/(2*l_k^2));
%Create the plot
figure(6)
plot(T,THETA);
%Format the plot
xlabel('time (sec)');
ylabel('angle (rad)');
ylim([0,1.2]);
title('the Turning Angle vs. Time for Explicit Simulation');
% toc