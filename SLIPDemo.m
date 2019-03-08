%% Spring Loaded Inverted Pendulumn (SLIP)
%  This is a non-disspative SLIP model 
%  Modified from the framework written by by Russ Tedrake.
%  ZC March 6,2019

%% Initialization
clear;clc;
% model properties
m   = 1;    % mass of body;
r0  = 1;    % initial length of leg
k   = 800;  % stiffness of spring
b   = 0;    % damping of spring was 8
g   = 9.8;  % gravity

dx_dsr = 2; % desired forward velocity ( m/s)
y_dsr = 1;% desired height ( meter )
K_raibert = 0.02; % raibert controller velocity gain
theta_max = pi/6;

% initiate model state
  x_td = 0;   % touchdown position
  angle_td = 0;  % touchdown angel
% in flight phase using Cartisian coordinate
  y0=y_dsr+dx_dsr^2/2/g;
  q0=[0;      % x horizontal position
     0;       % dx horizontal velocity
     y0;      % y vertical position
     0;       % dy vertical velocity
     1];      % phase flight = 1 || stance = 0

% time setting
t_start=0;    % start time
t_max=8;     % time at end of the simulation 
t_latest=t_start; % the beginning of each loop
t_step=1e-2;  % maximum step in ode45
t_td=0;       % touchdown time
t_lo=0;       % lift off time 
% output 
T = t_start;  % one column
Q = q0.';     % five columns 
P = q0(5)*angle_td+(1-q0(5))*x_td; % store angle or touchdown position 
% animation setting
flag=1;       % enable 1 || disable 0
flag_e=0;
%% Simulation
while(t_latest<t_max)
    
    tspan=[t_latest,t_max];
    
    
    if(q0(5))
        % flight phase = true 
        % set event function for ode45        
        options = odeset('Events',@(t,q) Event_Flight(t,q,r0,angle_td), ...
                    'MaxStep', t_step);
        % solve the dynamic equations
        [t,q,te,qe,ie] = ode45(@(t,q) Dynamics(q,r0,m,g,k,x_td,b), tspan, q0, options);
%         [t,q,te,qe,ie] = ode45(@(t,q) FlightDynamics(q,g), tspan, q0, options);
        % Record data
        T = [T; t];
        Q = [Q; q]; 
        p=ones(size(t,1),1).*(q0(5)*angle_td+(1-q0(5))*x_td);
        P = [P;p]; 
        % Update state
        if(size(ie,1))
            t_td=te+t_latest;
            t_latest = te;
            q0=qe;
            x_td=qe(1)+qe(3)*tan(angle_td); 
        %   q0=Flight2Stance(q_end,angle_td);
            q0(5)=0;
            disp('switch to stance  phase');
        else
            break;
        end
    else
        % stance phase 
        % set event function for ode45
         options = odeset('Events',@(t,q) Event_Stance(t,q,r0,x_td), ...
                    'MaxStep', t_step);
        % solve the dynamic equations
        [t,q,te,qe,ie] = ode45(@(t,q) Dynamics(q,r0,m,g,k,x_td,b), tspan, q0, options);
%         [t,q,te,qe,ie] = ode45(@(t,q) StanceDynamics(q,r0,m,g,k,x_td), tspan, q0, options);
        % Record data
        T = [T; t];
        Q = [Q; q]; 
        p=ones(size(t,1),1).*(q0(5)*angle_td+(1-q0(5))*x_td);
        P = [P;p];
        % Update state
        if(ie)
            t_lo=te+t_latest;
            t_latest = te;
            q0=qe;
    %         q0=Stance2Flight(q_end,x_td);
            % stance time
            t_stance = t_lo-t_td;
            % foot placement in next loop
            x_ft = q0(2)*t_stance/2 + K_raibert*(q0(2)-dx_dsr);
            % next touch down angle
            angle_td = asin(x_ft/r0);
            if angle_td > theta_max
                angle_td = theta_max;
            elseif angle_td < -theta_max
                angle_td = -theta_max;
            end
            q0(5)=1;
            disp('switch to flight  phase');
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
    
    if Q(j,5)
        r=r0;
    else
        r=sqrt(Q(j,3).^2+(Q(j,1)-P(j))^2);
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
    if(Q(i,5))
        toe=hip+r0*[sin(P(i)); -cos(P(i))];
    else
        toe=[P(i);0];
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
% plot(Q(:,3),Q(:,4));
% function xp=Flight2Stance(x,angle_td)
% % transform from Cartisian to Polar
% 
%     r = x(3)/cos(angle_td);
%     xp = [ r ; angle_td;
%       -x(2)*sin(angle_td) + x(4)*cos(angle_td);
%       -(x(2)*cos(angle_td) + x(4)*sin(angle_td))/r;];
% end
% 
% function x=Stance2Flight(xp,x_td)
%     x = [ x_td-xp(1)*sin(xp(3));...
%           xp(1)*cos(xp(3)); ...
%           -xp(2)*sin(xp(3)) - xp(1)*xp(4)*cos(xp(3)); ...
%           xp(2)*cos(xp(3)) - xp(1)*xp(4)*sin(xp(3))];
% end

% function dq=FlightDynamics(q,g)
%      % q = [x;xdot;y;ydot] in Cartisian coordinate
%     dq=[q(2);
%         0;
%         q(4);
%         -g];
% end

function dq=Dynamics(q,r0,m,g,k,x_td,b)
    if(q(5))
        r=r0;
    else
        r=sqrt((q(1)-x_td)^2+q(3)^2);
    end
    % q = [x;xdot;y;ydot] in Cartisian coordinate
    dq = [q(2);
          (q(1)-x_td)*k*((r0/r)-1)/m;       % CoM > touchdown position : accelerate
          q(4);
          q(3)*k*((r0/r)-1)/m-g;
          0];
end
% % function dq=StanceDynamics(q,r0,m,g,k)
% %      % q = [r;rdot;theta;theta_dot] in polar coordinate
% %      dq = [q(2);
% %            k*(r0-q(1))/m+q(1)*q(4)^2-g*cos(q(3));
% %            q(4);
% %            g*sin(q(3))/q(1)-2*q(2)*q(4)/q(1)];
% % end

function [zeroCrossing,isterminal,direction] = Event_Stance(t,q,r0,x_td)
    r=sqrt((q(1)-x_td)^2+q(3)^2);
    liftoff = r-r0;
    % Locate the time when height passes through zero in a 
    % decreasing direction and stop integration.
    zeroCrossing = liftoff; 
    isterminal   = 1;             
    direction    = 1;
end

function [zeroCrossing,isterminal,direction] = Event_Flight(t,q,r0,theta)
    touchdown    = r0*cos(theta) - q(3);
    % Locate the time when height passes through zero in a 
    % decreasing direction and stop integration.
    zeroCrossing = touchdown; 
    isterminal   = 1;             
    direction    = 1;
end