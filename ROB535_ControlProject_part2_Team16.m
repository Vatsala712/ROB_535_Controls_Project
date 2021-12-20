function [sol_2, Flag_terminate] = ROB535_ControlProject_part2_Team16(TestTrack,Xobs_seen,curr_state)
    % ROB535 Controls Project Team 16
    % Authors: Lars Grantz and Santoshi Kulkarni

    cline = TestTrack.cline;
    rBound = TestTrack.br;
    lBound = TestTrack.bl;

    min_track_width = min(vecnorm(lBound - rBound));
    iter = 2;
    input_range=[-0.5, 0.5;
                 -5000,5000];
    tracking_radius = min_track_width/2;
    x0 = curr_state;
    
    gains.kpd = 1.4;        %1.0    %1.4
    gains.kdd = 0.05;       %0.05   %0.05
    gains.kobsd = 0.77;      %0.4   %0.75
    gains.kpF = 11;         %10     %11
    gains.kobsF = 1.0;      %1.0    %1.0
    obs_radius = 7.5;        %10      %7
    max_iterations = 50;
    terminate_flag = false;
    dt = 0.01;
    
    cline_idx = closest_idx(x0, cline);
    
    % Getting inital Inputs
    u0 = get_u0(cline, cline_idx, x0, gains.kpd, input_range);
    %u0 = [0,0];
    x(1,:) = x0;
    u(1,:) = u0;
    u(2,:) = u0;
    x(2,:) = x(1,:) + dzdt_calc(x(1,:),u(2,:))'*dt;
    finish_count = 0;
    
    while(~terminate_flag && iter <= max_iterations)
        % update centerline goal point
        [cline_idx, curr_dist] = path_update(x,iter,cline,cline_idx, tracking_radius);

        % Calculating obstacle forces
        obstacle_Forces = get_obstacle_forces(Xobs_seen, x, iter, obs_radius);
        
        % Calculates next set of inputs
        u(iter+1,:) = get_next_inputs(cline, cline_idx, x, iter, gains, input_range, obstacle_Forces, curr_dist);

        % Update position using forward simulate
        x(iter+1,:) = x(iter,:) + dzdt_calc(x(iter,:),u(iter+1,:))'*dt;
        
        iter = iter + 1;
    end
    
    % Checks if needing to terminate
    T = 0:dt:(dt*size(x,1)-dt);
    info = getTrajectoryInfo(x,u,Xobs_seen,T,TestTrack);
    [terminate_flag] = check_terminate(info);
    
    percent_complete = info.percent_of_track_completed
    sol_2 = u;
    Flag_terminate = terminate_flag;
end

function cline_idx = closest_idx(x0, cline)
    % Finding closest cline index
    cline_idx = 1;
    curr_pos = [x0(1), x0(3)]';
    curr_dist = norm(curr_pos-cline(:,cline_idx));
    dist_to_iter = curr_dist;
    prev_dist = curr_dist + 1;
    
    width_idx = cline_idx;
    while(width_idx < size(cline,2) && prev_dist > dist_to_iter)
        width_idx = width_idx + 1;
        prev_dist = dist_to_iter;
        dist_to_iter = norm(curr_pos-cline(:,width_idx));
        if dist_to_iter > 50
            dist_to_iter = prev_dist-0.00001;
        end
    end
    cline_idx = width_idx;
end

function u0 = get_u0(cline, cline_idx, x0, kpd, input_range)
    % Gets the initial control inputs
    diff_x = cline(1,cline_idx) - x0(1);
    diff_y = cline(2,cline_idx) - x0(3);
    goal_angle = atan2(diff_y, diff_x);
    diff_theta = -x0(5) + goal_angle;
    
    delta0 = kpd*diff_theta;
    
    % Used dzdt equations for equilibrium assuming delta = 0
    F_x0 = 68.642 - 700*x0(4)*x0(6);
    
    % Check max inputs
    if abs(delta0) > input_range(1,2)
        delta0 = input_range(1,2)*sign(delta0);
    end
    if abs(F_x0) > input_range(2,2)
        F_x0 = input_range(2,2)*sign(F_x0);
    end
    u0 = [delta0,F_x0];
end

function [cline_idx, curr_dist] = path_update(x,iter,cline,cline_idx, tracking_radius)
    % This gets the next closest point along the centerline to follow
    curr_pos = [x(iter,1), x(iter,3)]';
    curr_dist = norm(curr_pos-cline(:,cline_idx));
    while (curr_dist <= tracking_radius && cline_idx ~= size(cline,2))
        cline_idx = cline_idx + 1;
        if cline_idx > size(cline,2)
            cline_idx = size(cline,2);
        end
        curr_dist = norm(curr_pos-cline(:,cline_idx));
    end
end

function obstacle_Forces = get_obstacle_forces(Xobs_seen, x, iter, obs_radius)
    % Goes through each of the seen obstacles and creates a small control
    % disturbance depending depending on the 
    obstacle_Forces = [0;0];
    for i = 1:size(Xobs_seen,2)
        obj_center = [mean(Xobs_seen{i}(:,1)), mean(Xobs_seen{i}(:,2))];
        diff_x_obs = obj_center(1) - x(iter,1);
        diff_y_obs = obj_center(2) - x(iter,3);
        dist_to_obs = norm([diff_x_obs, diff_y_obs]);
        obs_angle = atan2(diff_y_obs, diff_x_obs);
        diff_obs = -x(iter,5) + obs_angle;

        if dist_to_obs < obs_radius && abs(diff_obs) < pi/2
            obstacle_Forces(1) = obstacle_Forces(1) + 1*sign(diff_obs);
            obstacle_Forces(2) = obstacle_Forces(2) + 1/(dist_to_obs)^2;
        end
    end
end

function u_next = get_next_inputs(cline, cline_idx, x, iter, gains, input_range, obstacle_Forces, curr_dist)
    % Calculate Target Angle
    diff_x = cline(1,cline_idx) - x(iter,1);
    diff_y = cline(2,cline_idx) - x(iter,3);
    goal_angle = atan2(diff_y, diff_x);
    diff_theta = -x(iter,5) + goal_angle;

    % Apply & Store controls to the input
    diff_psi = x(iter,5) - x(iter-1,5);
    delta = gains.kpd*(diff_theta) - gains.kdd*diff_psi - gains.kobsd*(obstacle_Forces(1));
    F_x = gains.kpF*(curr_dist) - gains.kobsF*(obstacle_Forces(2));

    if cline_idx == size(cline,2)
        F_x = 100;
    end
    
    % Check Max Inputs
    if abs(delta) > input_range(1,2)
        delta = input_range(1,2)*sign(delta);
    end
    if abs(F_x) > input_range(2,2)
        F_x = input_range(2,2)*sign(F_x);
    end
    u_next = [delta, F_x];
end

function [terminate_flag] = check_terminate(info)
    terminate_flag = false;
    if (~isempty(info.left_track_position) || ~isempty(info.crash_position) || ~isempty(info.t_finished))
        terminate_flag = true;
    end
end

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