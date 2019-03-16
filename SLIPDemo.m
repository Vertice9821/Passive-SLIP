%% Spring Loaded Inverted Pendulumn (SLIP)
%  This is a non-disspative SLIP model 
%  Modified from the framework written by by Russ Tedrake.
%  ZC March 6,2019

%% Initialization
clear;clc;
% model properties
m   = 1;    % mass of body;
r0  = 0.8;    % initial length of leg
k   = 1500;  % stiffness of spring
g   = 9.8;  % gravity

dx_dsr = 1; % desired forward velocity ( m/s)
K_raibert = 0.02; % raibert controller velocity gain
theta_max = pi/8;

% initiate model state
% in flight phase using Cartisian coordinate
  q0=[0;      % x horizontal position
     0;     % dx horizontal velocity
     1;     % y vertical position
     0;       % dy vertical velocity
     0;       % phase flight = 0 || stance = 1
     0;       % x_td     touch down position 
     0];      % angle_td touch down angle from the vertical line

% time setting
t_start=0;    % start time
t_max=8;      % end of the simulation 
t_latest=t_start; % the beginning of each loop
t_step=1e-2;  % maximum step in ode45
t_td=0;       % touchdown time
t_lo=0;       % lift off time 
% output 
T = t_start;  % one column
Q = q0.';     % seven columns 
% animation setting
flag=0;       % enable 1 || disable 0
flag_e=~flag;
%% Simulation
% set event function for ode45
  options = odeset('Events',@(t,q) Event(t,q,r0), ...
                    'MaxStep', t_step);
                
while(t_latest<t_max)
    
    tspan=[t_latest,t_max];

    if(q0(5))
        % stance phase 
        % solve the dynamic equations
        [t,q,te,qe,ie] = ode45(@(t,q) Dynamics(q,r0,m,g,k), tspan, q0, options);
        % Record data
        T = [T; t];
        Q = [Q; q];
        % Update state at lift off
        if(ie)
            t_lo=te+t_latest;  % lift off time w.r.t. the whole time span
            t_latest = te;     % beginning of the next loop
            q0=qe;             
            % stance time
            t_stance = t_lo-t_td;
            % foot placement in next loop
            x_ft = q0(2)*t_stance/2 + K_raibert*(q0(2)-dx_dsr);
            % next touch down angle
            angle_td = asin(x_ft/r0);
             if abs(angle_td) > theta_max
                angle_td = sign(angle_td)*theta_max;
            end
            q0(7)=angle_td;             % record touch down angle 
            q0(6)=0;                    % clear touch down position
            q0(5)=0;                   
            disp('switch to flight  phase');
        else
            break;            
        end     
    else
        % flight phase 
        % solve the dynamic equations
        [t,q,te,qe,ie] = ode45(@(t,q) Dynamics(q,r0,m,g,k), tspan, q0, options);
        % Record data
        T = [T; t];
        Q = [Q; q]; 
        % Update state at touch down
        if(size(ie,1))
            t_td=te+t_latest;
            t_latest = te;
            q0=qe;
            x_td=qe(1)+qe(3)*tan(qe(7)); 
            q0(6)=x_td;        % record touch down position
            q0(7)=0;           % clear touch down angle
            q0(5)=1;
            disp('switch to stance  phase');
        else
            break;
        end
    end
        
end
%% Plot Energy
if flag_e
M=size(T,1);
Ek=[];
Ep=[];
for j=1:M
    Ek=[Ek;0.5*m*(Q(j,2).^2+Q(j,4).^2)];
    
    if ~Q(j,5)
        r=r0;
    else
        r=sqrt(Q(j,3).^2+(Q(j,1)-Q(j,6))^2);
    end
    Ep=[Ep;m*g*Q(j,3)+0.5*k*(r-r0)^2];
end
plot(T,Ek);
hold on;
plot(T,Ep);
hold on
plot(T,Ek+Ep);
hold on
legend('Ek','Ep','E');
end
%% Animation
if flag
V = VideoWriter('Passive SLIP.avi');
V.FrameRate = length(T)/t_max;
open(V);

f=figure(1);
subplot(1,2,1);
title('\fontsize{10}\fontname{Arial Black}SLIP Animation (2D)');

subplot(1,2,2);
h=animatedline;
title('\fontsize{10}\fontname{Arial Black}State ');
xlabel("y (m)");
ylabel("dydt (m/s)");
axis([0.5 1.5 -4 4]);
axis square;
hold on;

N=size(T,1);
for i=1:N
    hip=[Q(i,1);Q(i,3)];
    if(~Q(i,5))
        toe=hip+r0*[sin(Q(i,7)); -cos(Q(i,7))];
    else
        toe=[Q(i,6);0];
    end
    figure(f);
    subplot(1,2,1);
    cla;
    axis([-0.5 2 -0.5 2]);
    axis square;
    hold on;
    % leg
    line([hip(1);toe(1)], [hip(2); toe(2)],'Color',[0.459 0.059 0.427],'LineWidth',2);
    % hip
    t = 0:0.1:2*pi;
    line(hip(1)+0.15*sin(t),hip(2)+0.15*cos(t),'Color',[0 0 0]);
    fill(hip(1)+0.15*sin(t),hip(2)+0.15*cos(t),[ 0.867 0.639 0 ]);
    
    line([-10,10],[0,0]);

    subplot(1,2,2);
    addpoints(h,Q(i,3),Q(i,4));
    
    drawnow;
    
    frame(i) = getframe(f); 
    writeVideo(V,frame(i));
end
close(V);
end

function dq=Dynamics(q,r0,m,g,k)
    
    if(~q(5))  % in flight phase
        Fx=0;
        Fy=0;
    else       % in stance phase
        r=sqrt((q(1)-q(6))^2+q(3)^2);
        s=(q(1)-q(6))/r;       % sin(theta) 
        c=q(3)/r;             % cos(theta)  
        F_spring=k*(r0-r);
        Fx=F_spring*s;
        Fy=F_spring*c;
    end
    % q = [x;xdot;y;ydot;phase] in Cartisian coordinate
    dq = [q(2);
          Fx/m;       % CoM > touchdown position : accelerate
          q(4);
          Fy/m-g;
          0;
          0;
          0];
end

function [zeroCrossing,isterminal,direction] = Event(t,q,r0)
    if(q(5))
        r=sqrt((q(1)-q(6))^2+q(3)^2);
        liftoff = r-r0;
        zeroCrossing = liftoff; 
    else
        thouchdown=r0*cos(q(7))-q(3);
        zeroCrossing = thouchdown;
    end
    % Locate the time when height passes through zero in a 
    % decreasing direction and stop integration.
    isterminal   = 1;             
    direction    = 1;
end
