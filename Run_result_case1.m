close all;
clear;
addpath('utilities');
Paths = @Trajectory;
wp =        [0    0   0;
             1    1   1;
             2    0   2;
             3   -1   1;
             4    0   0]'; 
Paths([],[],wp);
Stabilizer = @PIDController;
[t, st, QP] = Simulator(Paths, Stabilizer);


%% Functions
%Trajectory function
function [ des_st ] = Trajectory(time1, st, wp)
persistent wp0 Time_Trajectory allCoeffs multFact time distance
if nargin > 2
    wp0 = wp';
    n = size(wp0, 1) - 1;
    multFact = 8;
    avgSpeed = 2; %m/sec
    d = wp(:,2:end) - wp(:,1:end-1);
    distance = 1/(avgSpeed) * sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    Time_Trajectory = [0, cumsum(distance)];
    A = zeros(multFact * n); 
    b = zeros(multFact * n, 3);
    numPts = size(wp0,1);
    b(1,:) = wp0(1,:);
    rowP = 1;
    for ipositionPt = 2:numPts-1
        row  = rowP + 1;
        rowP = row + 1;
        b(row,:)  = wp0(ipositionPt,:);
        b(rowP,:) = wp0(ipositionPt,:);
    end
    b(rowP+1,:) = wp0(end,:);    
    % Constraint values for A matrix
    positionCon = [1 0 0 0 0 0 0 0;
        1 1 1 1 1 1 1 1];
    drvtvCon = [0 1 2 3 4 5 6 7;
        0 0 2 6 12 20 30 42;
        0 0 0 6 24 60 120 210];    
    continuousCon = [0 1 2 3 4 5 6 7 0 -1 0 0 0 0 0 0;
        0 0 2 6 12 20 30 42 0 0 -2 0 0 0 0 0;
        0 0 0 6 24 60 120 210 0 0 0 -6 0 0 0 0;
        0 0 0 0 24 120 360 840 0 0 0 0 -24 0 0 0;
        0 0 0 0 0 120 720 2520 0 0 0 0 0 -120 0 0;
        0 0 0 0 0 0 720 5040 0 0 0 0 0 0 -720 0];
    % Waypoint constraints in A
    rowEnd = 0;
    for iposition = 1:n
        colStart = (iposition-1)*multFact + 1;
        colEnd   = colStart + multFact - 1;
        rowStart = rowEnd + 1;
        rowEnd   = rowStart + 1;
        A(rowStart:rowEnd,colStart:colEnd) = positionCon;
    end
    nxtRow = rowEnd;
    for iDrv = 1:3
        nxtRow = nxtRow + 1;
        A(nxtRow, iDrv+1) = 1;
    end
    
    startRow = (nxtRow + 1);
    endRow   = startRow + size(drvtvCon,1) - 1;
    startCol = (n-1) * multFact + 1;
    endCol   = size(A,2);
    A(startRow:endRow,startCol:endCol) = drvtvCon;
    
    colEnd = 0;
    for iCon = 1:n-1
        startRow = endRow + 1;
        endRow   = startRow + size(continuousCon,1) - 1;
        colStart = colEnd + 1;
        colEnd   = colStart + size(continuousCon,2) - 1;
        A(startRow:endRow,colStart:colEnd) = continuousCon;
        colEnd = multFact * iCon;
    end    
    allCoeffs = A\b;    
else
    if(time1 > Time_Trajectory(end))
        des_st.position = wp0(end,:)';
        des_st.velocity = zeros(3,1);
        des_st.acceleration = zeros(3,1);
    else
        timeindex = find(Time_Trajectory <= time1,1,'last');
        
        if(time1 == 0)
            des_st.position = wp0(1,:)';
            des_st.velocity = zeros(3,1);
            des_st.acceleration = zeros(3,1);
        else
            if(timeindex > 1)
                time1 = time1 - Time_Trajectory(timeindex);
            end
            p = time1/distance(timeindex);
            coeffs = allCoeffs(((timeindex-1)*multFact + 1): (timeindex*8), :);
            des_st.position = ([1, p, p^2, p^3, p^4, p^5, p^6, p^7] * coeffs)';
            des_st.velocity = ([0, 1, 2*p, 3*p^2, 4*p^3, 5*p^4, 6*p^5, 7*p^6] * coeffs)';
            des_st.acceleration = ([0, 0, 2, 6*p, 12*p^2, 20*p^3, 30*p^4, 42*p^5] * coeffs)';
        end
    end
    des_st.yw = 0;
    des_st.ywd = 0;
end
end
%Controller Function
function [Force, Moment] = PIDController(t, st, desired_st, parameters)
desiredrv = 0; 
desiredpv = 0; 
xaxis_kp = 100;
xaxis_kd = 14; 
yaxis_kp = 100;
yaxis_kd = 14;
zaxis_kp = 100; 
zaxis_kd = 14; 
phiKp = 100;
phiKd = 0.5;
thetKp = 100;
thetKd = 0.5;
psKp = 100;
psKd = 0.5;
%% control for acceleration 
x_acceleration = desired_st.acceleration(1) + xaxis_kd*(desired_st.velocity(1) - st.velocity(1)) + xaxis_kp*(desired_st.position(1) - st.position(1)); 
y_acceleration = desired_st.acceleration(2) + yaxis_kd*(desired_st.velocity(2) - st.velocity(2)) + yaxis_kp*(desired_st.position(2) - st.position(2));
z_acceleration = desired_st.acceleration(3) + zaxis_kd*(desired_st.velocity(3) - st.velocity(3)) + zaxis_kp*(desired_st.position(3) - st.position(3));
desiredph   = 1/parameters.gr * (x_acceleration * sin(desired_st.yw) - y_acceleration * cos(desired_st.yw));
desiredth = 1/parameters.gr * (x_acceleration * cos(desired_st.yw) + y_acceleration * sin(desired_st.yw));
%% Thrust 
Force = parameters.mass * (z_acceleration + parameters.gr);
%% Moment
phiu   = phiKp * (desiredph - st.rotat(1)) + phiKd * (desiredrv - st.omg(1));
thetau = thetKp * (desiredth - st.rotat(2)) + thetKd * (desiredpv - st.omg(2));
psiu   = psKp * (desired_st.yw - st.rotat(3)) + psKd * (desired_st.ywd - st.omg(3));
Moment = [phiu;
    thetau;
    psiu];
end

%Simulator function
function [t_out, s_out, QP] = Simulator(Path, Stabilizer)
addpath('utils');
real_time = true;
max_time = 50;
parameters = sys_parameters;
disp('Initializing figures...');
h_fig = figure;
h_3d = gca;
axis equal
grid on
view(48.8, 25.8);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]')
AUVcolors = lines(256);
set(gcf,'Renderer','OpenGL')
disp('Setting initial conditions...');
tstep    = 0.01;
cstep    = 0.02;
max_iter = max_time/cstep;
nstep    = cstep/tstep;
time     = 0;
err = []; 
% Start and stop position
des_start = Path(0, []);
des_stop  = Path(inf, []);
stop_position  = des_stop.position;
x0    = init_st(des_start.position, 0);
xtraj = zeros(max_iter*nstep, length(x0));
ttraj = zeros(max_iter*nstep, 1);
x       = x0;
position_tol = 0.01;
velocity_tol = 0.01;
disp('Simulation Running....');

for iter = 1:max_iter

    timeint = time:tstep:time+cstep;

    tic;

    % Initialize AUV plot
    if iter == 1
        QP = AUVPlot(1, x0, 0.1, 0.04, AUVcolors(1,:), max_iter, h_3d);
        current_st = stToQd(x);
        des_st = Path(time, current_st);
        QP.UpdateAUVPlot(x, [des_st.position; des_st.velocity], time);
        h_title = title(sprintf('iteration: %d, time: %4.2f', iter, time));
    end

    % Run simulation
    [tsave, xsave] = ode45(@(t,s) AUVEOM(t, s, Stabilizer, Path, parameters), timeint, x);
    x    = xsave(end, :)';

    % Save to traj
    xtraj((iter-1)*nstep+1:iter*nstep,:) = xsave(1:end-1,:);
    ttraj((iter-1)*nstep+1:iter*nstep) = tsave(1:end-1);

    % Update AUV plot
    current_st = stToQd(x);
    des_st = Path(time + cstep, current_st);
    QP.UpdateAUVPlot(x, [des_st.position; des_st.velocity], time + cstep);
    set(h_title, 'String', sprintf('iteration: %d, time: %4.2f', iter, time + cstep))

    time = time + cstep; % Update simulation time
    t = toc;
    % Check to make sure ode45 is not timing out
    if(t> cstep*50)
        err = 'Ode45 Unstable';
        break;
    end
    % Pause to make real-time
    if real_time && (t < cstep)
        pause(cstep - t);
    end
    % Check termination criteria
    if terminate_check(x, time, stop_position, position_tol, velocity_tol, max_time)
        break
    end
end
% Truncate xtraj and ttraj
xtraj = xtraj(1:iter*nstep,:);
ttraj = ttraj(1:iter*nstep);

% Truncate saved variables
QP.TruncateHist();
% Plot position
h_position = figure('Name', ['AUV position']);
plot_st(h_position, QP.st_hist(1:3,:), QP.time_hist, 'position', 'vic');
plot_st(h_position, QP.st_des_hist(1:3,:), QP.time_hist, 'position', 'des');
legend('Actual position','Desired position')

% Assuming QP.st_hist and QP.st_des_hist are matrices of size (3 x num_points)
% and QP.time_hist is a vector containing corresponding timestamps

% Calculate position error between the two sets of data
position_error = QP.st_hist(1:3,:) - QP.st_des_hist(1:3,:);

% Create a figure for the position error plot
h_position_error = figure('Name', 'position Error');

% Plot the x-axis position error
subplot(3, 1, 1);
plot(QP.time_hist, position_error(1,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('X-axis position Error');
grid on;

% Plot the y-axis position error
subplot(3, 1, 2);
plot(QP.time_hist, position_error(2,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('Y-axis position Error');
grid on;

% Plot the z-axis position error
subplot(3, 1, 3);
plot(QP.time_hist, position_error(3,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m]');
title('Z-axis position Error');
grid on;


% Plot velocity
h_velocity = figure('Name', ['AUV velocity']);
plot_st(h_velocity, QP.st_hist(4:6,:), QP.time_hist, 'velocity', 'vic');
plot_st(h_velocity, QP.st_des_hist(4:6,:), QP.time_hist, 'velocity', 'des');
legend('Actual velocity','Desired velocity')


% Assuming QP.st_hist and QP.st_des_hist are matrices of size (3 x num_points)
% and QP.time_hist is a vector containing corresponding timestamps

% Calculate velocity error between the two sets of data
velocity_error = QP.st_hist(4:6,:) - QP.st_des_hist(4:6,:);

% Create a figure for the velocity error plot
h_velocity_error = figure('Name', 'velocity Error');

% Plot the x-axis velocity error
subplot(3, 1, 1);
plot(QP.time_hist, velocity_error(1,:), 'r', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('X-axis velocity Error');
grid on;

% Plot the y-axis velocity error
subplot(3, 1, 2);
plot(QP.time_hist, velocity_error(2,:), 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('Y-axis velocity Error');
grid on;

% Plot the z-axis velocity error
subplot(3, 1, 3);
plot(QP.time_hist, velocity_error(3,:), 'b', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error [m/s]');
title('Z-axis velocity Error');
grid on;


if(~isempty(err))
end
disp('finished.')
t_out = ttraj;
s_out = xtraj;
end