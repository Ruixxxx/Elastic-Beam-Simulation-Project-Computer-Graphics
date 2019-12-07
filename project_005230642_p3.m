%%Generalized case of elastic beam falling in viscous flow
%Use the implicit Euler method to calculate the position and velocity of 
%N nodes of a beam in viscous fluid vs. time with initial conditions. Plot 
%the position and velocity of the middle node vs. time.

%%Typically scripts begin with
%{
Rui Xu
005230642
December 14,2018
%}
clc;clear all;close all;

%Assign the variables
N=21;           %number of nodes
r0=0.001;       %cross-sectional radius of the beam
E=1e9;          %Young's modulus
EA=E*pi*r0^2;   %stretching stiffness
EI=E*pi*r0^4/4; %bending stiffness
l=0.1;          %length of the beam
l_k=l/(N-1);    %length of each segment between two nodes

mu=1000;        %viscosity of the fluid
rho_fluid=1000; %density of the fluid
rho_metal=7000; %density of the sphere
R=ones(N,1)*l_k/10; %radii of the spheres
R((N+1)/2)=0.025;
[W1,m1]=WeightMass(R(1),rho_metal,rho_fluid); %weights and masses
[W_mid,m_mid]=WeightMass(R((N+1)/2),rho_metal,rho_fluid);
W=ones(N,1)*W1;
W((N+1)/2)=W_mid;
m=ones(N,1)*m1;
m((N+1)/2)=m_mid;

%%Implicit simulation
%Time counter
dt=1e-2;
tfinal=50;
%Calculate the steps in for-loop
tstep=ceil(tfinal/dt);

%Initialize the positions and velocities of the spheres
x=zeros(N,1);
y=zeros(N,1);
x_dot=zeros(N,1);
y_dot=zeros(N,1);

%Intialize the positions and velocities of the middle node
x_mid=zeros(tstep+1,1);
y_mid=zeros(tstep+1,1);
x_mid_dot=zeros(tstep+1,1);
y_mid_dot=zeros(tstep+1,1);

%Initial conditions
for k=1:N
    x(k)=l_k*(k-1);x_dot(k)=0;
    y(k)=0;y_dot(k)=0;
end

x_mid(1)=x((N+1)/2);
y_mid(1)=y((N+1)/2);
x_mid_dot(1)=x_dot((N+1)/2);
y_mid_dot(1)=y_dot((N+1)/2);

%Create the video
v=VideoWriter('BeamSimulationProblem3.avi');
open(v);

%Calculate the positions and velocities of the nodes implicitly
for k=1:tstep
    %Initialize the error
    error=10;
    %Initialize the tolerance
    tol=1e-3;
    %Set the starting points
    q=zeros(2*N,1);
    for i=1:N
        odd=2*i-1;
        even=2*i;
        q(odd)=x(i);
        q(even)=y(i);
    end
    while error>tol
        %Initialize the functions
        f=zeros(2*N,1);
        Fs=zeros(2*N,1);
        Fb=zeros(2*N,1);
        %Initialize the Jacobian
        J=zeros(2*N,2*N);
        Js=zeros(2*N,2*N);
        Jb=zeros(2*N,2*N);
        %Compute the functions and the Jacobian
        for i=1:N-1
            Fs_s=2*(i-1)+1;
            Fs_e=Fs_s+3;
            Fs(Fs_s:Fs_e,1)=Fs(Fs_s:Fs_e,1)+...
                gradEs(q(Fs_s),q(Fs_s+1),q(Fs_s+2),q(Fs_s+3),l_k,EA);
            Js(Fs_s:Fs_e,Fs_s:Fs_e)=Js(Fs_s:Fs_e,Fs_s:Fs_e)+...
                hessEs(q(Fs_s),q(Fs_s+1),q(Fs_s+2),q(Fs_s+3),l_k,EA)';
        end
        for i=2:N-1
            Fb_s=2*(i-1)-1;
            Fb_e=Fb_s+5;
            Fb(Fb_s:Fb_e,1)=Fb(Fb_s:Fb_e,1)+...
                gradEb(q(Fb_s),q(Fb_s+1),q(Fb_s+2),q(Fb_s+3),q(Fb_s+4),q(Fb_s+5),l_k,EI);
            Jb(Fb_s:Fb_e,Fb_s:Fb_e)=Jb(Fb_s:Fb_e,Fb_s:Fb_e)+...
                hessEb(q(Fb_s),q(Fb_s+1),q(Fb_s+2),q(Fb_s+3),q(Fb_s+4),q(Fb_s+5),l_k,EI)';
        end
        for i=1:N
            odd=2*i-1;
            even=2*i;
            f(odd)=m(i)/dt*((q(odd)-x(i))/dt-x_dot(i))+...
                (6*pi*mu*R(i))*((q(odd)-x(i))/dt);
            f(even)=m(i)/dt*((q(even)-y(i))/dt-y_dot(i))+...
                W(i)+(6*pi*mu*R(i))*((q(even)-y(i))/dt);
            J(odd,odd)=m(i)/(dt^2)+(6*pi*mu*R(i))/dt;
            J(even,even)=m(i)/(dt^2)+(6*pi*mu*R(i))/dt;
        end
        f=f+Fs+Fb;
        J=J+Js+Jb;
        %Newton's update
        dq=J\f;
        q=q-dq;
        %Compute the error
        error=sum(abs(f));
    end
    %Velocity
    for i=1:N
        odd=2*(i-1)+1;
        even=2*i;
        x_dot(i)=(q(odd)-x(i))/dt;
        y_dot(i)=(q(even)-y(i))/dt;
    end
    %Position
    for i=1:N
        odd=2*(i-1)+1;
        even=2*i;
        x(i)=q(odd);
        y(i)=q(even);
    end

    %Plot the shape of beam
    figure(1)
    plot(x,y,'.-');
    xlabel('x (m)');
    ylabel('y (m)');
    title(num2str(k*dt, 'time %4.2f s'));
    axis equal
    drawnow
    
    %Save frame to video
    writeVideo(v,getframe(gcf));
    
    x_mid(k+1)=x((N+1)/2);
    y_mid(k+1)=y((N+1)/2);
    x_mid_dot(k+1)=x_dot((N+1)/2);
    y_mid_dot(k+1)=y_dot((N+1)/2);
end

%Finish the video shooting
close(v);

%%Plot the positions and velocities of the middle node vs. time
t=0:dt:tfinal;
%Create the plot
figure(2)
plot(t,x_mid);
hold on
plot(t,y_mid);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('position (m)');
legend('x-position','y-position','Location','Best');
title('the Middle Node Position for Implicit Simulation');
%Create the plot
figure(3)
plot(t,x_mid_dot);
hold on
plot(t,y_mid_dot);
hold off
%Format the plot
xlabel('time (sec)');
ylabel('velocity (m/sec)');
legend('x-velocity','y-velocity','Location','Best');
title('the Middle Node Velocity for Implicit Simulation');

%%Display the terminal velocities
fprintf('The terminal velocity along x-axis is %f.\n',x_mid_dot(tstep+1));
fprintf('The terminal velocity along y-axis is %f.\n',y_mid_dot(tstep+1));

% %%Plot the terminal velocity vs. dt and N
% dt=[1,0.5,0.1,0.05,0.01,0.005,0.001];
% N=[3,5,11,21,51];
% TerminalVelocity=[-0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833;
%                  -0.005927,-0.005857,-0.005837,-0.005834,-0.005833];
% surf(N,dt,TerminalVelocity);
% ay=gca;
% ay.YScale='log';
% xlabel('number of nodes');
% ylabel('timestep size');
% zlabel('terminal velocity');
% title('The Terminal Velocity vs. Number of Nodes and Timestep Size');