%%Elastic beam bending
%Use the implicit Euler method to simulate the beam vs. time with initial
%conditions and force and plot the maximum vertical displacement vs. time.

%%Typically scripts begin with
%{
Rui Xu
005230642
December 14,2018
%}
clc;clear all;close all;

%%Assign the variables
N=50;        %number of nodes
l=1;         %length of the beam
l_k=l/(N-1); %length of each segment between two nodes
R=0.013;     %outer radius
r=0.011;     %inner radius
P=-2000;     %force
d=0.75;      %distance between the location of the applied force and the left-hand edge
E=7e10;      %modulus of elasticity
I=pi/4*(R^4-r^4); %moment of inertia of the cross section
EI=E*I;      %bending stiffness
EA=E*pi*(R^2-r^2); %stretching stiffness
rho=2700;    %the density of aluminum
m=pi*(R^2-r^2)*l*rho/(N-1); %mass of nodes

%%Implicit simulation
%Time counter
dt=1e-2;
tfinal=1;
%Calculate the steps in for-loop
tstep=ceil(tfinal/dt);

%Initialize the positions of the nodes
x=zeros(N,1);
y=zeros(N,1);
x_dot=zeros(N,1);
y_dot=zeros(N,1);

%Intialize the maximum vertical displacement
y_max=zeros(tstep+1,1);

%Initial conditions
for k=1:N
    x(k)=l_k*(k-1);x_dot(k)=0;
    y(k)=0;y_dot(k)=0;
end

y_max(1)=0;

%Create the video
v=VideoWriter('BeamSimulationProblem4.avi');
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
        %Find the location of the force
        [distance,Ind]=min(abs(x-d));
        f=Fs+Fb;
        f(1)=q(1);
        f(2)=q(2);
        f(2*N)=q(2*N);
        f(Ind*2)=f(Ind*2)-P;
        J=Js+Jb;
        J(1,:)=0;
        J(1,1)=1;
        J(2,:)=0;
        J(2,2)=1;
        J(2*N,:)=0;
        J(2*N,2*N)=1;
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
    
    %Find the maximum vertical displacement
    y_max(k+1)=(-1)*max(abs(y));
end

%Finish the video shooting
close(v);

%%Plot the maximum vertical displacement vs. time
t=0:dt:tfinal;
%Create the plot
figure(2)
plot(t,y_max);
%Format the plot
xlabel('time (sec)');
ylabel('y (m)');
title('the Maximum Vertical Displacement of the Beam');

% %%Calculate the theoretical solution
% c=min(d,l-d);
% ymax_theo=P*c*(l^2-c^2)^1.5/(9*sqrt(3)*E*I*l);
% %%Display the maximum displacements
% fprintf('The theoretical maximum displacement is %f.\n',ymax_theo);
% fprintf('The maximum displacement from implicit simulation is %f.\n',y_max(tstep+1));
% 
% %%Plot the maximum vertical displacements vs. P
% P=[-20000,-15000,-10000,-5000,-2000,-1000,-500];
% y_max=[-0.199924,-0.182167,-0.148229,-0.088935,-0.037109,-0.018674,-0.009352];
% ymax_theo=[-0.380449,-0.285337,-0.190225,-0.095112,-0.038045,-0.019022,-0.009511];
% figure(3)
% plot(P,y_max);
% hold on
% plot(P,ymax_theo);
% hold off
% xlabel('force (N)');
% ylabel('maximum displacement (m)');
% legend('implicit simulation','beam theory','Location','Best');
% title('the Maximum Vertical Displacement vs. Force');
