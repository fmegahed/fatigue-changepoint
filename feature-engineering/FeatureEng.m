clear all; close all; clc


%%
SamplingRate = 51.2;
%% data read, process and save segments
% Read data
imu_read = importdata('ShimmerData.xlsx');
imu_data = imu_read.data;
Data_in = imu_data(1200:1500,2:7); 

StepLengthAvg = 60;
PeakAccAvg = 5;

% Data processing and segmentation
[M_i_k_shimmer, Seg_points_S, Seg_points_E] = process_segment_3h(Data_in,SamplingRate,StepLengthAvg,PeakAccAvg);
save('MIK_SpecDef','M_i_k_shimmer')

shimmerSteps = M_i_k_shimmer(:,1:2);

figure('Color','w') 
lenSave = 0;
for i = 1:size(shimmerSteps,1) 
    
    plot(shimmerSteps{i,1}+lenSave+0.01,shimmerSteps{i,2}, 'Color',[0.2,0.2,0.6] ,'LineWidth',2)
    
    lenSave = lenSave + range(shimmerSteps{i,1});
    hold on
    stepStats(i,3) = range(shimmerSteps{i,1});
    stepStats(i,4) = range(shimmerSteps{i,2});
end
%axis equal
set(gca, 'FontSize', 20)
title('Transformed Gait Step Position Trajectories Derived from the Velocity Data','FontSize', 24)
xlabel('X (m)')
ylabel('Y (m)')
legend({'{\bf{P}}_{Frenet-Serret}'},'Location','southwest')


hold on
    
%%

function [M_i_k, Seg_points_S, Seg_points_E] = process_segment_3h(Data_in,SamplingRate,StepLengthAvg,PeakAccAvg)

    %% Read Data
    % This code imports the down sampled acceleration and gyro data and uses the butterworth low pass filter 
    % to smooth the transformed data from Kalman filter and segment the data into gait steps 
    acc_data = Data_in(:,1:3);
    gyro_data = Data_in(:,4:6);
    duration = (length(Data_in)-1)/SamplingRate;
    dt = 1/SamplingRate;
    t = (0:dt:duration)';
    %% Kalman filter 
    [State, Rot] = IMU_Kalman(Data_in, SamplingRate); 
    filtered_angle = State(:,7:9);      % filtered angle
    unbiased_gyro = State(:,10:12);     % unbiased gyro 
    
    t_plot = linspace(0,size(acc_data,1)/SamplingRate,size(acc_data,1));

    

    % finds the quaternions of the rotation in order to transfering
    % body frame set of data to global frame
    acc_glob = Quaternion(acc_data,Rot);
    gyro_glob = Quaternion(gyro_data,Rot);
    
    
    acc_g = repmat([0 -9.55 0],length(acc_data),1);
    acc_lin = acc_glob + acc_g; 
    %% Segmentation and Kinematics Calculation
    % Butterworth filter parameters for acceleration 
    n=4;fc=4;Fs=100;
    Wn = (2/Fs)*fc;
    [b,a]=butter(n,Wn,'low');
    Acceleration_filt = filter(b,a,acc_lin);
    Acceleration_filt_magn = filter(b,a,sqrt(sum(Acceleration_filt.^2,2)));
    
    
    % Segmentation based on peaks in acc data
    Seg_points_S = [];
    Seg_points_E = [];
    
    n = 21;
    while n<length(t)-StepLengthAvg
        max_b = max(Acceleration_filt_magn(n+fix(StepLengthAvg/2):n+StepLengthAvg));
        max_a = max(Acceleration_filt_magn(n:n+fix(StepLengthAvg/2)));
        if (max_b/max_a>1.2 || max_a/max_b>1.2) && max(Acceleration_filt_magn(n:n+StepLengthAvg))>PeakAccAvg      
            [~,I1] = max(Acceleration_filt_magn(n:n+StepLengthAvg));
            if n+I1+20 <= length(Acceleration_filt_magn)  
                [~,I2] = min(Acceleration_filt_magn(n+I1:n+I1+20));
                [~,I3] = min(Acceleration_filt_magn(n-20:n+10));
                if abs(diff([n-20+I3, n+I1+I2])-StepLengthAvg) < 0.3*StepLengthAvg
                    Seg_points_S = [Seg_points_S; n-20+I3];
                    Seg_points_E = [Seg_points_E; n+I1+I2-2];
                    n = n + I1+I2;
                else
                    n = n + 1;
                end
            else
                n = n + 1;
            end
        else
            n = n + 1;
        end
    end
    
    % Making segments windows
    k = 1;
    for j = 1:length(Seg_points_E)
        S_s(k) = Seg_points_S(j);
        S_e(k) = Seg_points_E(j);
        k = k + 1;
    end     
    
    figure('Color','w') 
    plot(t_plot,Acceleration_filt_magn,'Color',[0.2,0.2,0.6], 'LineWidth',2)
    hold on 
    plot(t_plot(S_s), Acceleration_filt_magn(S_s), 'o', 'MarkerSize',10,'MarkerEdgeColor',[.1 .6 .6],'MarkerFaceColor',[.1 .6 .6])
    hold on
    plot(t_plot(S_e), Acceleration_filt_magn(S_e), '^', 'MarkerSize',10,'MarkerEdgeColor',[1 .7 .6],'MarkerFaceColor',[1 .7 .6])
    set(gca, 'FontSize', 20)
    title('Acceleration Magnitude Data and Gait Segments Post Applying the Butterworth Filter','FontSize', 24)
    xlabel('Time (sec)')
    ylabel('Acceleration (m/s^2)')
    legend({'{\bf{a}}_{magnitude}','Gait Cycle Start Point','Gait Cycle Start Point'},'Location','northwest')

    axs = [1 2 3];

    for direc = axs

        % Velocity calculation
        Vel{1,direc} = zeros(length(t),1);
        for i = 1:length(S_s)-1
            Res_acc = trapz(t(S_s(i):S_e(i)),acc_lin(S_s(i):S_e(i),direc))/range(t(S_s(i):S_e(i))); 
            Vel{1,direc}(S_s(i):S_e(i),:) = cumtrapz(t(S_s(i):S_e(i)),acc_lin(S_s(i):S_e(i),direc)) - Res_acc*linspace(0,t(S_e(i))-t(S_s(i)),length(t(S_s(i):S_e(i))))';
        end

        % Position calculation
        Pos{1,direc} = cumtrapz(t,Vel{1,direc});
% 
        
        % Jerk calculation
        for i = 1:length(t)-1
            Jrk{1,direc}(i,:) = diff(acc_lin(i:i+1,direc))/dt;
        end

    end

    % Metrics Segmentation
    vel_xyz = cell2mat(Vel);
    Velocity = mag(vel_xyz,1);
    Jerk = mag(cell2mat(Jrk),1);
    
    
    figure('Color','w')
    subplot(3,1,1)
    plot(t_plot(1:100),acc_data(1:100,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(1:100),acc_data(1:100,2), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(1:100),acc_data(1:100,3), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
    set(gca, 'FontSize', 20)
    ylabel([{'IMU Acceleration'}; {'(m/s^2)'}])
    legend({'\bf{a_x}','\bf{a_x}','\bf{a_x}'},'Location','southwest')
    set(gca,'xtick',[])

    subplot(3,1,2)
    plot(t_plot(1:100),gyro_data(1:100,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(1:100),gyro_data(1:100,2), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(1:100),gyro_data(1:100,3), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
    set(gca, 'FontSize', 20)
    ylabel([{'IMU Angular'}; {'Velocity (deg/s)'}])
    legend({'\theta_{\bf{x}}','\theta_{\bf{y}}','\theta_{\bf{z}}'},'Location','southwest')
    set(gca,'xtick',[])
    

    subplot(3,1,3)
    plot(t_plot(20:100),vel_xyz(20:100,1), t_plot(100:120),zeros(21,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(20:100),vel_xyz(20:100,2), t_plot(100:120),zeros(21,1), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
    hold on
    plot(t_plot(20:100),vel_xyz(20:100,3), t_plot(100:120),zeros(21,1), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
    set(gca, 'FontSize', 20)
    xlabel('Time (sec)')
    ylabel([{'Undrifted Foot'}, {'Velocity (m/s)'}])
    legend({'\bf{V_x}','\bf{V_x}','\bf{V_x}'},'Location','southwest')
    
    figure('Color','w')
    plot(t_plot(20:45),zeros(26,1), t_plot(45:100),...
        Acceleration_filt_magn(45:100)-linspace(Acceleration_filt_magn(45),Acceleration_filt_magn(100),...
        length(Acceleration_filt_magn(45:100)))', t_plot(100:120),zeros(21,1), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
    set(gca, 'FontSize', 20)
    ylabel([{'Foot Acceleration'}, {'Magnitude (m/s^2)'}])
    legend('{\bf{a}}_{magnitude}','Location','southwest')
    set(gca,'xtick',[])
    
    m = [Pos{1,1}(2:end), Pos{1,2}(2:end), Pos{1,3}(2:end), t(2:end), Velocity(2:end), t(2:end), Acceleration_filt_magn(2:end), t(2:end), Jerk, filtered_angle(2:end,1), filtered_angle(2:end,2), filtered_angle(2:end,1), unbiased_gyro(2:end,1), filtered_angle(2:end,2), unbiased_gyro(2:end,2),filtered_angle(2:end,3), unbiased_gyro(2:end,3), t(2:end)];
    
    
    SegStep = Seg_points_E-Seg_points_S;
    
    for i = 1:size(m,2)
        k = 1;
        for j = 1:length(SegStep)
            % Making the data start from zero in each segment window
            if i == size(m,2) 
                M_i_k{k,i} = m(S_s(j):S_e(j),i);
            else
                M_i_k{k,i} = m(S_s(j):S_e(j),i) - m(S_s(j),i);
            end
            k = k + 1;
        end
    end
    
   
    for i = 1:length(M_i_k(:,1))
        x = M_i_k{i,1};
        y = M_i_k{i,2};
        z = M_i_k{i,3};
        
        sf = fit([x, y],z,'poly11');
        
        x_SP = ([x(end),y(end),z(end)]-[x(1),y(1),z(1)])./norm([x(end),y(end),z(end)]-[x(1),y(1),z(1)]);
        z_SP = [sf.p10,sf.p01,-1]./norm([sf.p10,sf.p01,-1]);  
        y_SP = cross(x_SP,z_SP);

        M = [x_SP;y_SP;z_SP];
        xyz_SP = M*[x,y,z]';
        
        M_i_k{i,1} = xyz_SP(1,:);
        M_i_k{i,2} = xyz_SP(2,:);
    end
    
    M_i_k(:,3)=[];    
    
    
end

function [State, R] = IMU_Kalman(imu_data, SamplingRate)

duration = (length(imu_data)-1)/SamplingRate;
dt = 1/SamplingRate;
t = (0:dt:duration)';

acc_data = imu_data(:,1:3);
gyro_data = imu_data(:,4:6);
gyro = gyro_data;

v_gyro = sqrt(0.3)*randn(length(t),3); % measurement noise for gyro,  variance = 0.3
v_acc = sqrt(0.3)*randn(length(t),3); % measurement noise for accellerometer, variance = 0.3

% Compute the angles computed by using only accelerometers of gyroscope
x = [-atan(acc_data(:,1)./sqrt(sum([acc_data(:,2).^2 acc_data(:,3).^2],2)))*(180/pi), ...
     atan(acc_data(:,2)./sqrt(sum([acc_data(:,1).^2 acc_data(:,3).^2],2)))*(180/pi), ...
     atan(acc_data(:,3)./sqrt(sum([acc_data(:,1).^2 acc_data(:,2).^2],2)))*(180/pi)];
x_acc = x + v_acc; % Angle computed from accelerometer measurement

t_plot = linspace(0,size(x,1)/SamplingRate,size(x,1));

figure('Color','w') 
plot(t_plot,x(:,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
hold on
plot(t_plot,x(:,2), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
hold on
plot(t_plot,x(:,3), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
set(gca, 'FontSize', 20)
title('Angular Position Calculated from Accelerometer','FontSize', 24)
xlabel('Time (sec)')
ylabel('Angle (degree)')
legend({'{\theta}_x','{\theta}_y','{\theta}_z'},'Location','southeast')


x_gyro = cumsum(gyro*dt,1); % Angle computed by integration of gyro measurement

figure('Color','w')
plot(t_plot,x_gyro(:,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
hold on
plot(t_plot,x_gyro(:,2), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
hold on
plot(t_plot,x_gyro(:,3), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
set(gca, 'FontSize', 20)
title('Angular Position Calculated from Gyro','FontSize', 24)
xlabel('Time (sec)')
ylabel('Angle (degree)')
legend({'{\theta}_x','{\theta}_y','{\theta}_z'},'Location','southeast')



P = [0.1 0; 0 0.1];  
R_angle = 0.5;   

Q_angle = 0.9;
Q_gyro = 0.1;
Q = [Q_angle 0; 0 Q_gyro];

A = [0 1; 0 1];
q_bias = [0 0 0]; % Initialize gyro bias
angle = [-4 81 -8]; % Initialize gyro angle
q_m = 0;
X = [-4 81 -8; 0 0 0];

    for i=1:length(t)

         % Gyro update 

         q_m = gyro(i,:);

         q = q_m - q_bias; % gyro bias removal

         Pdot = A*P + P*A' + Q;

         rate = q;

         angle = angle + q*dt;

         P = P + Pdot*dt;

         % Kalman (Accelerometer) update 

         H = [0 1]; 
         angle_err = x_acc(i,:)-angle;
         E = H*P*H' + R_angle;

         K = P*H'*inv(E);

         P = P - K*H*P;
         X = X + K * angle_err;
         
         
         x1(i,:) = X(1,:);
         R(:,:,i) = ang2orth(x1(i,:));
         x2(i,:) = X(2,:);
         angle = x1(i,:);
         q_bias = x2(i,:);

         x3(i,:) = q;  % unbiased gyro rate
    end

State = [x_acc x_gyro x1 x3];

x_acc = State(:,1:3);   % angle calculated from accelerometer
x_gyro = State(:,4:6);  % angle calculated from gyro 
filtered_angle = State(:,7:9);      % filtered angle
unbiased_gyro = State(:,10:12);     % unbiased gyro
rotation_mat = R;


figure('Color','w') ;
plot(t_plot,filtered_angle(:,1), '--','Color',[1,0.4,0.6], 'LineWidth',2')
hold on
plot(t_plot,filtered_angle(:,2), '-.','Color',[0,0.9,0.6], 'LineWidth',2')
hold on
plot(t_plot,filtered_angle(:,3), '-', 'Color',[0.2,0.8,1], 'LineWidth',2')
set(gca, 'FontSize', 20)
title('Angular Position Results from the Kalman Filter','FontSize', 24)
xlabel('Time (sec)')
ylabel('Angle (degree)')
legend({'{\theta}_x','{\theta}_y','{\theta}_z'},'Location','southeast')


    function orthm = ang2orth(ang) 
        sa = sind(ang(1)); ca = cosd(ang(1)); 
        sb = sind(ang(2)-80); cb = cosd(ang(2)-80); 
        sc = sind(ang(3)); cc = cosd(ang(3)); 

        ra = [  ca,  sa,  0; ... 
               -sa,  ca,  0; ... 
                 0,   0,  1]; 
        rb = [  cb,  0, -sb; ... 
                 0,  1,  0; ... 
                sb,  0,  cb]; 
        rc = [  1,   0,   0; ... 
                0,   cc, sc;... 
                0,  -sc, cc]; 
        orthm = (ra*rb*rc)';
    end

end



function X_glob = Quaternion(X_body,R)

Q = dcm2q(R);
X_glob = qvqc(Q,X_body);

%% Functions


    function q=dcm2q(R)
        % DCM2Q(R) converts direction cosine matrices into quaternions.
        %
        %     The resultant quaternion(s) will perform the equivalent vector
        %     transformation as the input DCM(s), i.e.:
        %
        %       qconj(Q)*V*Q = R*V
        %
        %     where R is the DCM, V is a vector, and Q is the quaternion.  Note that
        %     for purposes of quaternion-vector multiplication, a vector is treated
        %     as a quaterion with a scalar element of zero.
        %
        %     If the input is a 3x3xN array, the output will be a vector of
        %     quaternions where input direction cosine matrix R(:,:,k) corresponds
        %     to the output quaternion Q(k,:).
        %
        %     Note that this function is meaningless for non-orthonormal matrices!
        %
        % See also Q2DCM.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.11 $
        % $Date: 2009-07-25 04:28:18 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.

        % Thanks to Tatsuki Kashitani for pointing out the numerical instability in
        % the original implementation.  His suggested fix also included a check for
        % the "sr" values below being zero.  But I think I've convinced myself that
        % this isn't necessary if the input matrices are orthonormal (or at least
        % very close to orthonormal).

        if nargin~=1
          error('One input argument required');
        else
          size_R=size(R);
          if ( size_R(1)~=3 || size_R(2)~=3 || length(size_R)>3 )
            error('Invalid input: must be a 3x3xN array')
          end
        end

        q = zeros( 4, size( R, 3 ) );

        for id_dcm = 1 : size( R, 3 )
          dcm = R( :, :, id_dcm );
          if trace( dcm ) > 0
            % Positve Trace Algorithm
            sr  = sqrt( 1 + trace( dcm ));
            sr2 = 2*sr;
            q(1,id_dcm) = ( dcm(2,3) - dcm(3,2) ) / sr2;
            q(2,id_dcm) = ( dcm(3,1) - dcm(1,3) ) / sr2;
            q(3,id_dcm) = ( dcm(1,2) - dcm(2,1) ) / sr2;
            q(4,id_dcm) = 0.5 * sr;
          else
            % Negative Trace Algorithm
            if ( dcm(1,1) > dcm(2,2) ) && ( dcm(1,1) > dcm(3,3) )
              % Maximum Value at DCM(1,1)
              sr  = sqrt( 1 + (dcm(1,1) - ( dcm(2,2) + dcm(3,3) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = 0.5 * sr;
              q(2,id_dcm) = ( dcm(2,1) + dcm(1,2) ) / sr2;
              q(3,id_dcm) = ( dcm(3,1) + dcm(1,3) ) / sr2;
              q(4,id_dcm) = ( dcm(2,3) - dcm(3,2) ) / sr2;
            elseif dcm(2,2) > dcm(3,3)
              % Maximum Value at DCM(2,2)
              sr  = sqrt( 1 + (dcm(2,2) - ( dcm(3,3) + dcm(1,1) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = ( dcm(2,1) + dcm(1,2) ) / sr2;
              q(2,id_dcm) = 0.5 * sr;
              q(3,id_dcm) = ( dcm(2,3) + dcm(3,2) ) / sr2;
              q(4,id_dcm) = ( dcm(3,1) - dcm(1,3) ) / sr2;
            else
              % Maximum Value at DCM(3,3)
              sr  = sqrt( 1 + (dcm(3,3) - ( dcm(1,1) + dcm(2,2) )) );
              sr2 = 2*sr;
              q(1,id_dcm) = ( dcm(3,1) + dcm(1,3) ) / sr2;
              q(2,id_dcm) = ( dcm(2,3) + dcm(3,2) ) / sr2;
              q(3,id_dcm) = 0.5 * sr;
              q(4,id_dcm) = ( dcm(1,2) - dcm(2,1) ) / sr2;
            end
          end
        end

        % Make quaternion vector a column of quaternions
        q=q.';

        q=real(q);
    end


    function qtype=isq(q)
        % ISQ(Q) checks to see if Q is a quaternion or set of quaternions.
        %     ISQ returns a value accordingly:
        %
        %        0 if Q is not a quaternion or vector of quaternions:
        %          has more than 2 dimensions or neither dimension is of length 4
        %       
        %        1 if the component quaternions of Q are column vectors:
        %          Q is 4xN, where N~=4, or
        %          Q is 4x4 and only the columns are normalized 
        %
        %        2 if the component quaternions of Q are row vectors:
        %          Q is Nx4, where N~=4, or
        %          Q is 4x4 and only the rows are normalized 
        %
        %        3 if the shape of the component quaternions is indeterminant:
        %          Q is 4x4, and either both the columns and rows are normalized
        %          or neither the columns nor rows are normalized.
        %
        %     In other words, if Q is 4x4, ISQ attempts to discern the shape of
        %     component quaternions by determining whether the rows or the columns
        %     are normalized (i.e., it assumes that normalized quaternions are
        %     the more typical use of quaternions).
        %
        %     The test for normalization uses 2*EPS as a tolerance.
        %
        % See also ISNORMQ, EPS.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.7 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.

        if nargin~=1

          error('isq() requires one input argument');

        else

          tol=2*eps;

          size_q=size(q);

          if ( length(size_q)~=2 || max(size_q==4)~=1 )
            qtype=0; % Not a quaternion or quaternion vector

          elseif ( size_q(1)==4 && ...
                   ( size_q(2)~=4 || ( ~sum((sum(q.^2,1)-ones(1,4))>tol) &&   ...
                                        sum((sum(q.^2,2)-ones(4,1))>tol)    ) ...
                     ) ...
                   )
            qtype=1; % Component q's are column vectors

          elseif ( size_q(2)==4 && ...
                   ( size_q(1)~=4 || ( ~sum((sum(q.^2,2)-ones(4,1))>tol) &&   ...
                                        sum((sum(q.^2,1)-ones(1,4))>tol)    ) ...
                     ) ...
                   )
            qtype=2; % Component q's are row vectors

          else
            qtype=3; % Component q's are either columns or rows (indeterminate)

          end

        end
    end


    function qout=qconj(qin)
        % QCONJ(Q) calculates the conjugate of the quaternion Q.
        %     Works on "vectors" of quaterions as well.  Will return the same shape
        %     vector as input.  If input is a vector of four quaternions, QCONJ will
        %     determine whether the quaternions are row or column vectors according
        %     to ISQ.
        %
        % See also ISQ.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.16 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=1
          error('qconj() requires one input argument');
        else
          qtype = isq(qin);
          if ( qtype==0 )
            error(...
              'Invalid input: must be a quaternion or a vector of quaternions')
          elseif ( qtype==3 )
            warning(...
              'qconj:indeterminateShape', ...
              'Component quaternion shape indeterminate, assuming row vectors')
          end
        end

        % Make sure component quaternions are row vectors
        if( qtype == 1 )
          qin=qin.';
        end

        qout(:,1)=-qin(:,1);
        qout(:,2)=-qin(:,2);
        qout(:,3)=-qin(:,3);
        qout(:,4)= qin(:,4);

        % Make sure output is same shape as input
        if( qtype == 1 )
          qout=qout.';
        end
    end

    function v_out=qcvq(q,v)
        % QcVQ(Q,V) performs the operation qconj(Q)*V*Q
        %     where the vector is treated as a quaternion with a scalar element of
        %     zero.
        %
        %     Q and V can be vectors of quaternions and vectors, but they must
        %     either be the same length or one of them must have a length of one.
        %     The output will have the same shape as V.  Q will be passed through
        %     QNORM to ensure it is normalized.
        %
        % See also QVQc, QNORM, QMULT.

        % Note that QNORM is invoked by QMULT, therefore QcQV does not invoke
        % it directly.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.2 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('Two input arguments required.');
        else

          qtype=isq(q);
          if ( qtype == 0 )
            error('Input Q must be a quaternion or a vector of quaternions')
          end

          size_v=size(v);
          if ( length(size_v)~=2 || max(size_v==3)~=1 )
            error(['Invalid input: second input must be a 3-element vector', 10, ...
                   'or a vector of 3-element vectors'])
          end

        end

        % Make sure q is a column of quaternions
        if ( qtype==1 )
          q=q.';
        end

        % Make sure v is a column of vectors
        row_of_vectors = (size_v(2) ~= 3);
        if ( row_of_vectors )
          v=v.';
          size_v=size(v);
        end

        size_q=size(q);

        if (  size_q(1)~=size_v(1) && size_q(1)~=1 && size_v(1)~=1 )
          error(['Inputs do not have the same number of elements:', 10, ...
                 '   number of quaternions in q = ', num2str(size_q(1)), 10,...
                 '   number of vectors in v     = ', num2str(size_v(1)), 10,...
                 'Inputs must have the same number of elements, or', 10, ...
                 'one of the inputs must have a single element.']) 
        elseif ( size_q(1)==1 && size_v(1)==3 )
          if ( qtype==1 )
            warning(...
              'qcvq:assumingVcols', ...
              'Q is 4x1 and V is 3x3: assuming vectors are column vectors')
            row_of_vectors = 1;
            v=v.';
          else
            warning(...
              'qcvq:assumingVrows', ...
              'Q is 1x4 and V is 3x3: assuming vectors are row vectors')
          end
        elseif ( qtype==3 && size_v(1)==1 )
          if ( row_of_vectors )
            warning(...
              'qcvq:assumingQcols', ...
              'Q is 4x4 and V is 3x1: assuming quaternions are column vectors')
            q=q.';
          else
            warning(...
              'qcvq:assumingQrows', ...
              'Q is 4x4 and V is 1x3: assuming quaternions are row vectors')
          end  
        end

        % Build up full vectors if one input is a singleton
        if ( size_q(1) ~= size_v(1) )
          ones_length = ones(max(size_q(1),size_v(1)),1);
          if ( size_q(1) == 1 )
            q = [q(1)*ones_length ...
                 q(2)*ones_length ...
                 q(3)*ones_length ...
                 q(4)*ones_length ];
          else % size_v(1) == 1
            v = [v(1)*ones_length ...
                 v(2)*ones_length ...
                 v(3)*ones_length ];    
          end
        end

        % Add an element to V
        v(:,4)=zeros(size_v(1),1);

        % Turn off warnings before calling qconj (it has simillar warnings as
        % qvxform, so all warnings would just be duplicated).  Save current state of
        % warnings, though.
        warning_state = warning; warning('off', 'qconj:indeterminateShape');
        local_warning = lastwarn;

        % Perform transform
        vt=qmult(qconj(q),qmult(v,q));

        % Restore warning state to original state
        warning(warning_state);
        lastwarn(local_warning);

        % Eliminate last element of vt for output
        v_out = vt(:,1:3);

        % Make sure output vectors are the same shape as input vectors
        if ( row_of_vectors )
          v_out = v_out.';
        end
    end
    
    
    function q_out=qmult(q1,q2)
        % QMULT(Q1,Q2) calculates the product of two quaternions Q1 and Q2.
        %    Inputs can be vectors of quaternions, but they must either have the
        %    same number of component quaternions, or one input must be a single
        %    quaternion.  QMULT will determine whether the component quaternions of
        %    the inputs are row or column vectors according to ISQ.
        %  
        %    The output will have the same shape as Q1.  If the component
        %    quaternions of either Q1 or Q2 (but not both) are of indeterminate
        %    shape (see ISQ), then the shapes will be assumed to be the same for
        %    both inputs.  If both Q1 and Q2 are of indeterminate shape, then both
        %    are assumed to be composed of row vector quaternions.
        %
        % See also ISQ.

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.14 $
        % $Date: 2009-07-26 20:05:12 $

        % Copyright (c) 2001-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('qmult() requires two input arguments');
        else
          q1type = isq(q1);
          if ( q1type == 0 )
            error(['Invalid input: q1 must be a quaternion or a vector of' ...
                  ' quaternions'])
          end
          q2type = isq(q2);
          if ( q2type == 0 )
            error(['Invalid input: q2 must be a quaternion or a vector of' ...
                  ' quaternions'])
          end
        end

        % Make sure q1 is a column of quaternions (components are rows)
        if ( q1type==1 || (q1type==3 && q2type==1) )
          q1=q1.';
        end

        % Make sure q2 is a column of quaternions (components are rows)
        if ( q2type==1 || (q2type==3 && q1type==1) )
          q2=q2.';
        end

        num_q1=size(q1,1);
        num_q2=size(q2,1);

        if (  num_q1~=num_q2 && num_q1~=1 && num_q2~=1 )
          error(['Inputs do not have the same number of elements:', 10, ...
                 '   number of quaternions in q1 = ', num2str(num_q1), 10,...
                 '   number of quaternions in q2 = ', num2str(num_q2), 10,...
                 'Inputs must have the same number of elements, or', 10, ...
                 'one of the inputs must be a single quaternion (not a', 10, ...
                 'vector of quaternions).']) 
        end

        % Build up full quaternion vector if one input is a single quaternion
        if ( num_q1 ~= num_q2 )
          ones_length = ones(max(num_q1,num_q2),1);
          if ( num_q1 == 1 )
            q1 = [q1(1)*ones_length ...
                  q1(2)*ones_length ...
                  q1(3)*ones_length ...
                  q1(4)*ones_length ];
          else % num_q2 == 1
            q2 = [q2(1)*ones_length ...
                  q2(2)*ones_length ...
                  q2(3)*ones_length ...
                  q2(4)*ones_length ];    
          end
        end

        % Products

        % If q1 and q2 are not vectors of quaternions, then:
        %
        %   q1*q2 = q1*[ q2(4) -q2(3)  q2(2) -q2(1)
        %                q2(3)  q2(4) -q2(1) -q2(2)
        %               -q2(2)  q2(1)  q2(4) -q2(3)
        %                q2(1)  q2(2)  q2(3)  q2(4) ]
        %
        % But to deal with vectorized quaternions, we have to use the ugly
        % commands below.

        prod1 = ...
            [ q1(:,1).*q2(:,4) -q1(:,1).*q2(:,3)  q1(:,1).*q2(:,2) -q1(:,1).*q2(:,1)];
        prod2 = ...
            [ q1(:,2).*q2(:,3)  q1(:,2).*q2(:,4) -q1(:,2).*q2(:,1) -q1(:,2).*q2(:,2)];
        prod3 = ...
            [-q1(:,3).*q2(:,2)  q1(:,3).*q2(:,1)  q1(:,3).*q2(:,4) -q1(:,3).*q2(:,3)];
        prod4 = ...
            [ q1(:,4).*q2(:,1)  q1(:,4).*q2(:,2)  q1(:,4).*q2(:,3)  q1(:,4).*q2(:,4)];

        q_out = prod1 + prod2 + prod3 + prod4;

        % Make sure output is same format as q1
        if ( q1type==1 || (q1type==3 && q2type==1) )
          q_out=q_out.';
        end

        % NOTE that the following algorithm proved to be slower than the one used
        % above:
        %
        % q_out = zeros(size(q1));
        % 
        % q_out(:,1:3) = ...
        %     [q1(:,4) q1(:,4) q1(:,4)].*q2(:,1:3) + ...
        %     [q2(:,4) q2(:,4) q2(:,4)].*q1(:,1:3) + ...
        %     cross(q1(:,1:3), q2(:,1:3));
        % 
        % q_out(:,4) = q1(:,4).*q2(:,4) - dot(q1(:,1:3), q2(:,1:3), 2);
    end

    
    function v_out=qvqc(q,v)
        % QVQc(Q,V) performs the operation Q*V*qconj(Q)
        %     where the vector is treated as a quaternion with a scalar element of
        %     zero. 
        %
        %     Q and V can be vectors of quaternions and vectors, but they must
        %     either be the same length or one of them must have a length of one.
        %     The output will have the same shape as V.  Q will be passed through
        %     QNORM to ensure it is normalized.
        %
        % See also QcQV, QNORM

        % Release: $Name: quaternions-1_3 $
        % $Revision: 1.1 $
        % $Date: 2009-07-24 19:14:44 $

        % Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.


        if nargin~=2
          error('Two input arguments required');
        else
          q     = qconj(q);
          v_out = qcvq(q, v);
        end
    end
    
end


function N = mag(T,n)
% MAGNATUDE OF A VECTOR (Nx3)
%  M = mag(U)
N = sum(abs(T).^2,2).^(1/2);
d = find(N==0); 
N(d) = eps*ones(size(d));
N = N(:,ones(n,1));  
end

%}