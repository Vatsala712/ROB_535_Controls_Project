close all;
clear;
clc

load("TestTrack.mat")
figure('Position',[50 200 500 500]);
plot(TestTrack.cline(1,:), TestTrack.cline(2,:),'b:') %Length of 246
hold on
plot(TestTrack.bl(1,:), TestTrack.bl(2,:),'k')
plot(TestTrack.br(1,:), TestTrack.br(2,:),'k')
axis square
f = waitbar(0,'Starting','Name','Trials Completion');

time_avg = 0;
n = 1;  % Number of trials to run
complete_count = 0;
fail_count = 0;
j = 1;
dt = 0.01;
while j <= n %&& fail_count == 0
    [Y,U,t_total,t_update, Xobs] = forwardIntegrateTester();

    info = getTrajectoryInfo(Y,U,Xobs,(dt*size(Y,1)-dt),TestTrack)
    
    if ~isempty(info.t_finished)
        time_avg = time_avg + info.t_finished;
        complete_count = complete_count+1;
        complete_time(complete_count) = info.t_finished;
    else
        fail_count = fail_count+1;
        fail.x{fail_count} = Y;
        fail.u{fail_count} = U;
        fail.Xobs{fail_count} = Xobs;
        fail_percent(fail_count) = info.percent_of_track_completed;
    end
    
    waitbar(double(j)/n,f,sprintf('%3.0f Successes Rate With %1.0f Trials',double(complete_count)/j*100,j));
    j = j + 1;
end
close(f);

Success_Count = complete_count
if complete_count ~= 0
    time_avg = time_avg/complete_count
end

obj_list = zeros(size(Xobs,2),2);
for i = 1:size(Xobs,2)
    obj_list((4*(i-1)+1):(4*i),:) = Xobs{i};
end

plot(Y(:,1),Y(:,3),'r')
scatter(obj_list(:,1),obj_list(:,2),'k.')
if ~isempty(info.left_track_position)
    plot(info.left_track_position(1),info.left_track_position(2),'cx');
end
if ~isempty(info.crash_position)
    plot(info.crash_position(1),info.crash_position(2),'cx');
end
hold off;