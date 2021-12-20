close all;
clear;
clc;

load("TestTrack.mat")
figure(1);
plot(TestTrack.cline(1,:), TestTrack.cline(2,:),'b:') %Length of 246
hold on
plot(TestTrack.bl(1,:), TestTrack.bl(2,:),'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:),'k')
hold off;

% Loading TestTrack Data into variables
cline = TestTrack.cline;
rBound = TestTrack.br;
lBound = TestTrack.bl;
orientData = TestTrack.theta;
min_track_width = min(vecnorm(lBound - rBound));
nstates=6;
ninputs=2;

% Loop initialization / constants
x0 = [287,5,-176,0,2,0];
dt = 0.01;
input_range=[-0.5, 0.5;
             -5000,5000];
tracking_radius = min_track_width/2;
finish_radius = min_track_width/2;
terminate_flag = false;
finish_count = 0;
max_iterations = 25000;
iter = 2;
cline_idx = 1;
dt = 0.01;
Xobs = {};

F_x0 = 68.642 - 700*x0(4)*x0(6);
if abs(F_x0) > input_range(2,2)
    F_x0 = input_range(2,2)*sign(F_x0);
end
u0 = [0,0];
x(1,:) = x0;
u(1,:) = u0;
u(2,:) = u0;
x(2,:) = x(1,:) + dzdt_calc(x(1,:),u(2,:))'*dt;
kpd = 0.3; %0.3
kid = 0.0;
kdd = 0.01;
kobsd = 0.0;
kpF = 7; %10
kiF = 0.0;
kobsF = 0.0;
kid_sum = 0;
kiF_sum = 0;
f = waitbar(0,'Percent Completion');

tic
while(~terminate_flag && iter < max_iterations)
    % update centerline goal point
    curr_pos = [x(iter,1), x(iter,3)]';
    curr_dist = norm(curr_pos-cline(:,cline_idx));
    while (curr_dist <= tracking_radius && cline_idx ~= size(cline,2))
        cline_idx = cline_idx + 2;
        if cline_idx > size(cline,2)
            cline_idx = size(cline,2);
        end
        curr_dist = norm(curr_pos-cline(:,cline_idx));
    end
    
    % Calculate Target Angle
    diff_x = cline(1,cline_idx) - x(iter,1);
    diff_y = cline(2,cline_idx) - x(iter,3);
    goal_angle = atan2(diff_y, diff_x);
    diff_theta = -x(iter,5) + goal_angle;
    
    % Apply & Store controls to the input
    diff_psi = x(iter,5) - x(iter-1,5);
    delta = kpd*(diff_theta) + kid*(kid_sum) - kdd*diff_psi - kobsd*(0);
    kid_sum = kid_sum + diff_theta;
    F_x = kpF*(curr_dist) + kiF*(kiF_sum) - kobsF*(0);
    kiF_sum = kiF_sum + F_x;
    
    % Check Max Inputs
    if abs(delta) > input_range(1,2)
        delta = input_range(1,2)*sign(delta);
    end
    if cline_idx == size(cline,2)
        F_x = F_x * 10; 
    end
    if abs(F_x) > input_range(2,2)
        F_x = input_range(2,2)*sign(F_x);
    end
    u(iter+1,:) = [delta, F_x];
    %u(iter+1,:) = u(iter+1,:) - (u(iter+1,:) - u(iter,:))*0.75;
    
    % Update position using forward simulate
    %[x_new,T] = forwardIntegrateControlInput(u(iter:(iter+1),:),x(iter,:));
    x(iter+1,:) = x(iter,:) + dzdt_calc(x(iter,:),u(iter+1,:))'*dt;
    
    info = getTrajectoryInfo(x(iter:(iter+1),:),u(iter:(iter+1),:),Xobs,0:dt:0.01,TestTrack);
    if (~isempty(info.left_track_position) || ~isempty(info.t_finished))
        terminate_flag = true;
    end
        
    if mod(iter,100) == 0
        waitbar(info.percent_of_track_completed,f);
    end
    iter = iter + 1;
end
close(f);
toc
TrackTime = iter*dt;
iter
[x_check,T] = forwardIntegrateControlInput(u, x0);
info = getTrajectoryInfo(x_check,u,Xobs,T,TestTrack)
TrackTime = info.t_end;
left_track_position = info.left_track_position;
crash_position = info.crash_position;
percent_of_track_completed = info.percent_of_track_completed;

hold on;
plot(x(:,1),x(:,3),'r')
plot(x_check(:,1),x_check(:,3),'g--')
if ~isempty(info.left_track_position)
    plot(left_track_position(1),left_track_position(2),'cx');
end
if ~isempty(info.crash_position)
    plot(crash_position(1),crash_position(2),'cx');
end
hold off;

% Error = vecnorm(x_check'-x');
% figure(2);
% plot(T,Error)


function dzdt=dzdt_calc(x,U)
%constants
Nw=2;
f=0.01;
Iz=2667;
a=1.35;
b=1.45;
By=0.27;
Cy=1.2;
Dy=0.7;
Ey=-1.6;
Shy=0;
Svy=0;
m=1400;
g=9.806;


%generate input functions
delta_f=U(1);
F_x=U(2);

%slip angle functions in degrees
a_f=rad2deg(delta_f-atan2(x(4)+a*x(6),x(2)));
a_r=rad2deg(-atan2((x(4)-b*x(6)),x(2)));

%Nonlinear Tire Dynamics
phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

F_zf=b/(a+b)*m*g;
F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

F_zr=a/(a+b)*m*g;
F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

F_total=sqrt((Nw*F_x)^2+(F_yr^2));
F_max=0.7*m*g;

if F_total>F_max
    
    F_x=F_max/F_total*F_x;
  
    F_yr=F_max/F_total*F_yr;
end

%vehicle dynamics
dzdt= [x(2)*cos(x(5))-x(4)*sin(x(5));...
          (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+x(4)*x(6);...
          x(2)*sin(x(5))+x(4)*cos(x(5));...
          (F_yf*cos(delta_f)+F_yr)/m-x(2)*x(6);...
          x(6);...
          (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
end


