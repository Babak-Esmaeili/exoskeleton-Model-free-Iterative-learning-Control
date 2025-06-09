clc;
clear;
close all;

tic;
%% Time sequence
T0 = 0;
Ts = 0.0001;
Tf = 1;

t = T0:Ts:Tf;

nSteps = length(t);

%% Controller parameters
% MFAILITSMC
eta = 1;
mu = 1;

epsilon = 1e-5;

rho = 3/5;

c_1 = 10*diag([1,1,1]);
c_2 = 6*diag([0.1,0.1,0.1]);

c_s = 0.00001*diag([0.002,0.002,0.0002]);

m_1 = 0*diag([1,1,1]);
m_2 = 0.01*diag([0.1,1,20]);
m_3 = 0*[1 1 1]';
m = [m_1 m_2 m_3];

nIter = 50;

%%
% hip angle parameters
ah1 = 22.47;
bh1 = 5.467;
ch1 = 2.226;
ah2 = 5.492;
bh2 = 11.15;
ch2 = -1.265;
ah3 = 7.951; 
bh3 = 1.212;
ch3 = 1.939;
ah4 = 1.586;
bh4 = 19.16;
ch4 = -0.2865;

% knee angle parameters
ak1 = 119.9;
bk1 = 2.324;
ck1 = -1.488;
ak2 = 14.75;
bk2 = 12.16;
ck2 = 2.142;
ak3 = 108.1;
bk3 = 3.446;
ck3 = 1.23;
ak4 = 3.773;
bk4 = 18.5;
ck4 = 3.165;

% ankle angle parameters
aL1 = 15;
bL1 = 0.9133;
cL1 = 2.271;
aL2 = 8.752;
bL2 = 11.12;
cL2 = 3.413;
aL3 = 3.969;
bL3 = 20.06;
cL3 = -2.352;
aL4 = 1.219;
bL4 = 29.21;
cL4 = 4.468;

deg2rad = pi/180;

%% Initialization
nInputs = 3;
nOutputs = 3;

u = zeros(nInputs,nSteps,nIter);
u(:,:,1) = zeros(nInputs,nSteps);
delta_u_eq = zeros(nInputs,nSteps,nIter);
delta_u_sw = zeros(nInputs,nSteps,nIter);
delta_u = zeros(nInputs,nSteps,nIter);

y = zeros(nOutputs,nSteps+1,nIter);
y(:,:,1) = zeros(nOutputs,nSteps+1);
delta_y = zeros(nOutputs,nSteps+1,nIter);

delta_omega = zeros(nOutputs+nInputs+1,nSteps+1,nIter);

x = zeros(nOutputs,nSteps+1,nIter);

hip_angle = ah1*sin(bh1*t(1)+ch1) + ah2*sin(bh2*t(1)+ch2) + ah3*sin(bh3*t(1)+ch3) + ah4*sin(bh4*t(1)+ch4);
hip_angle_d = ah1*bh1*cos(bh1*t(1)+ch1) + ah2*bh2*cos(bh2*t(1)+ch2) + ah3*bh3*cos(bh3*t(1)+ch3) + ah4*bh4*cos(bh4*t(1)+ch4);
hip_angle_dd = -ah1*bh1*bh1*sin(bh1*t(1)+ch1) - ah2*bh2*bh2*sin(bh2*t(1)+ch2) - ah3*bh3*bh3*sin(bh3*t(1)+ch3) + ah4*bh4*bh4*sin(bh4*t(1)+ch4);

knee_angle = ak1*sin(bk1*t(1)+ck1) + ak2*sin(bk2*t(1)+ck2) + ak3*sin(bk3*t(1)+ck3) + ak4*sin(bk4*t(1)+ck4);
knee_angle_d = ak1*bk1*cos(bk1*t(1)+ck1) + ak2*bk2*cos(bk2*t(1)+ck2) + ak3*bk3*cos(bk3*t(1)+ck3) + ak4*bk4*cos(bk4*t(1)+ck4);
knee_angle_dd = -ak1*bk1*bk1*sin(bk1*t(1)+ck1) - ak2*bk2*bk2*sin(bk2*t(1)+ck2) - ak3*bk3*bk3*sin(bk3*t(1)+ck3) - ak4*bk4*bk4*sin(bk4*t(1)+ck4);

ankle_angle = aL1*sin(bL1*t(1)+cL1) + aL2*sin(bL2*t(1)+cL2) + aL3*sin(bL3*t(1)+cL3) + aL4*sin(bL4*t(1)+cL4);
ankle_angle_d = aL1*bL1*cos(bL1*t(1)+cL1) + aL2*bL2*cos(bL2*t(1)+cL2) + aL3*bL3*cos(bL3*t(1)+cL3) + aL4*bL4*cos(bL4*t(1)+cL4);
ankle_angle_dd = -aL1*bL1*bL1*sin(bL1*t(1)+cL1) - aL2*bL2*bL2*sin(bL2*t(1)+cL2) - aL3*bL3*bL3*sin(bL3*t(1)+cL3) - aL4*bL4*bL4*sin(bL4*t(1)+cL4);

y_1_d = hip_angle*deg2rad;
y_2_d = knee_angle*deg2rad;
y_3_d = ankle_angle*deg2rad;

y_1_dot_d = hip_angle_d*deg2rad;
y_2_dot_d = knee_angle_d*deg2rad;
y_3_dot_d = ankle_angle_d*deg2rad;

y_1_ddot_d = hip_angle_dd*deg2rad;
y_2_ddot_d = knee_angle_dd*deg2rad;
y_3_ddot_d = ankle_angle_dd*deg2rad;

y_d = zeros(nOutputs,nSteps+1);
y_d(:,1) = [y_1_dot_d y_2_dot_d y_3_dot_d]';

x_d = zeros(nOutputs,nSteps+1);

e = zeros(nOutputs,nSteps+1,nIter);

int_e = zeros(nOutputs,nSteps+1,nIter);

s = zeros(nOutputs,nSteps+1,nIter);
delta_s = zeros(nOutputs,nSteps+1,nIter);

PHI_hat = cell(nSteps,nIter);
PHI_hat_1 = cell(nSteps,nIter);
PHI_hat_2 = cell(nSteps,nIter);
PHI_hat_3 = cell(nSteps,nIter);

dist = zeros(nOutputs,nSteps);

e_max = zeros(nOutputs,nIter);
e_rms = zeros(nOutputs,nIter);
e_rms(:,1) = 3*[1 1 1]';
e_mean = zeros(nOutputs,nIter);
e_mean(:,1) = [0.006 0.002 0.002]';
e_IAE = zeros(nOutputs,nIter);
e_IAE(:,1) = [0.2056 0.1969 0.3760]';

fprintf('Simulating MFAILITSMC on 2-DOF exoskeleton...\n');
fprintf('Number of iterations is %d\n',nIter);
fprintf('-----------------------------------------------\n');

%% Iteration loop
for i = 2:nIter
    
    fprintf('Iteration #%d is evaluated...\n',i);
    
    q = [0 0 0]';
%     q_dot = 0.01*rand(3,1);
%     q_dot = [-3;0.01*rand(2,1)];
    q_dot = -0.01*rand(3,1);
%     q_dot = [0 0 0]';
    y(:,1,i) = q_dot;

    for k = 1:nSteps
        
        %% PJM estimation
        if(i == 2)
            
            PHI_hat_1{k,1} = 1.1*diag([1,1,1]);
            PHI_hat_2{k,1} = 0.5*diag([.005,.001,.01]);
            PHI_hat_3{k,1} = 1*[0 0 0]';
            
            PHI_hat{k,1} = [PHI_hat_1{k,1} PHI_hat_2{k,1} PHI_hat_3{k,1}];
            
        end
        
        PHI_hat{k,i} = PHI_hat{k,i-1} + eta*((delta_y(:,k+1,i-1)-(PHI_hat{k,i-1}+m)*delta_omega(:,k,i-1))*delta_omega(:,k,i-1)')/(mu+norm(delta_omega(:,k,i-1))^2);
        
        PHI_hat_1{k,i} = PHI_hat{k,i}(:,1:nOutputs);
        PHI_hat_2{k,i} = PHI_hat{k,i}(:,nOutputs+1:2*nOutputs);
        PHI_hat_3{k,i} = PHI_hat{k,i}(:,2*nOutputs+1);
        
        for ii = 1:nInputs
            for jj = 1:nOutputs
                if(sign(PHI_hat_2{k,i}(ii,jj)) ~= sign(PHI_hat_2{k,1}(ii,jj)) || norm(PHI_hat_2{k,i}) <= epsilon)    % better!
                    PHI_hat_2{k,i}(ii,jj) = PHI_hat_2{k,1}(ii,jj);
                end
            end
        end
        
        PHI_hat{k,i} = [PHI_hat_1{k,i} PHI_hat_2{k,i} PHI_hat_3{k,i}];
                
        %% Reference signals
        hip_angle = ah1*sin(bh1*t(k)+ch1) + ah2*sin(bh2*t(k)+ch2) + ah3*sin(bh3*t(k)+ch3) + ah4*sin(bh4*t(k)+ch4);
        hip_angle_d = ah1*bh1*cos(bh1*t(k)+ch1) + ah2*bh2*cos(bh2*t(k)+ch2) + ah3*bh3*cos(bh3*t(k)+ch3) + ah4*bh4*cos(bh4*t(k)+ch4);
        hip_angle_dd = -ah1*bh1*bh1*sin(bh1*t(k)+ch1) - ah2*bh2*bh2*sin(bh2*t(k)+ch2) - ah3*bh3*bh3*sin(bh3*t(k)+ch3) + ah4*bh4*bh4*sin(bh4*t(k)+ch4);

        knee_angle = ak1*sin(bk1*t(k)+ck1) + ak2*sin(bk2*t(k)+ck2) + ak3*sin(bk3*t(k)+ck3) + ak4*sin(bk4*t(k)+ck4);
        knee_angle_d = ak1*bk1*cos(bk1*t(k)+ck1) + ak2*bk2*cos(bk2*t(k)+ck2) + ak3*bk3*cos(bk3*t(k)+ck3) + ak4*bk4*cos(bk4*t(k)+ck4);
        knee_angle_dd = -ak1*bk1*bk1*sin(bk1*t(k)+ck1) - ak2*bk2*bk2*sin(bk2*t(k)+ck2) - ak3*bk3*bk3*sin(bk3*t(k)+ck3) - ak4*bk4*bk4*sin(bk4*t(k)+ck4);

        ankle_angle = aL1*sin(bL1*t(k)+cL1) + aL2*sin(bL2*t(k)+cL2) + aL3*sin(bL3*t(k)+cL3) + aL4*sin(bL4*t(k)+cL4);
        ankle_angle_d = aL1*bL1*cos(bL1*t(k)+cL1) + aL2*bL2*cos(bL2*t(k)+cL2) + aL3*bL3*cos(bL3*t(k)+cL3) + aL4*bL4*cos(bL4*t(k)+cL4);
        ankle_angle_dd = -aL1*bL1*bL1*sin(bL1*t(k)+cL1) - aL2*bL2*bL2*sin(bL2*t(k)+cL2) - aL3*bL3*bL3*sin(bL3*t(k)+cL3) - aL4*bL4*bL4*sin(bL4*t(k)+cL4);

        y_1_d = hip_angle*deg2rad;
        y_2_d = knee_angle*deg2rad;
        y_3_d = ankle_angle*deg2rad;

        y_1_dot_d = hip_angle_d*deg2rad;
        y_2_dot_d = knee_angle_d*deg2rad;
        y_3_dot_d = ankle_angle_d*deg2rad;

        y_1_ddot_d = hip_angle_dd*deg2rad;
        y_2_ddot_d = knee_angle_dd*deg2rad;
        y_3_ddot_d = ankle_angle_dd*deg2rad;

        x_d(:,k+1) = [y_1_d y_2_d y_3_d]';
        y_d(:,k+1) = [y_1_dot_d y_2_dot_d y_3_dot_d]';
    
        %% MFAILITSMC
        delta_u_eq(:,k,i) = ((PHI_hat_2{k,i}+m_2)^-1)*(y_d(:,k+1)-y(:,k+1,i-1)- ...
                            (PHI_hat_1{k,i}+m_1)*delta_y(:,k,i)-PHI_hat_3{k,i}+...
                            (c_1^-1)*c_2*int_e(:,k+1,i-1)-(c_1^-1)*s(:,k+1,i-1));
        
        delta_u_sw(:,k,i) = c_s*((PHI_hat_2{k,i}+m_2)^-1)*sign(s(:,k,i));
%         delta_u_sw(:,k,i) = c_s*((PHI_hat_2{k,i}+m_2)^-1)*tanh(s(:,k,i))/0.1;

        delta_u(:,k,i) = delta_u_eq(:,k,i) + delta_u_sw(:,k,i);
        u(:,k,i) = u(:,k,i-1) + delta_u(:,k,i);

        %% Plant simulation
        q_ddot = threeDOF_exoskeleton_dynamics(q,q_dot,u(:,k,i),dist(:,k));
        q_dot = q_dot + Ts*q_ddot;
        q = q + Ts*q_dot;
        
        x(:,k+1,i) = q;
        
        y(:,k+1,i) = q_dot;
        
        %% Tracking errors
        e(:,k+1,i) = y_d(:,k+1) - y(:,k+1,i);
        
        %% Iterative sliding surfaces
        s(:,k+1,i) = c_1*e(:,k+1,i) + c_2*int_e(:,k+1,i-1);
        
%         int_e(:,k+1,i) = int_e(:,k,i) + sig_func(e(:,k+1,i),rho);
        
        %% Difference of iterative sliding surfaces
        delta_s(:,k+1,i) = s(:,k+1,i) - s(:,k+1,i-1);
        
        %% Incremental data
        delta_y(:,k+1,i) = y(:,k+1,i) - y(:,k+1,i-1);
        
        delta_omega(:,k,i) = [delta_y(:,k,i);delta_u(:,k,i);1];
        
    end
    int_e(:,:,i) = int_e(:,:,i-1) + sig_func(e(:,:,i),rho);
    
    e_rms(:,i) = rms(e(:,:,i)');
    e_mean(:,i) = sum(e(:,:,i),2)/(nSteps+1);
    e_IAE(:,i) = sum(abs(e(:,:,i)),2)*Ts;
    e_max(:,i) =  max(abs(e(:,5:end,i)),[],2);
    
end

% s_last = zeros(nOutputs,nIter);
% my_step = nSteps;
% for ii = 1:nIter
% 
%     s_last(:,ii) = s(:,my_step,ii);
%     
% end

y_d = y_d(:,1:end-1);

y = y(:,1:end-1,:);

e = e(:,1:end-1,:);

s = s(:,1:end-1,:);
delta_s = delta_s(:,1:end-1,:);

% e_IAE(:,1) = (e_IAE(:,2)-e_IAE(:,3)) + e_IAE(:,2);  % good!

fprintf('\nSuccessfully done!\n');
fprintf('-----------------------------------------------\n');

toc;
%% Plot results
LineWidth = 1.5;
LineWidth_d = 2;

my_iter = 5;

%% Difference of iterative sliding surfaces in last iteration number
figure(1);
subplot(3,1,1);
plot(t,delta_s(1,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\Delta{s_1}');

title('Difference of sliding surfaces in last iteration number');

subplot(3,1,2);
plot(t,delta_s(2,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\Delta{s_2}');

subplot(3,1,3);
plot(t,delta_s(3,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\Delta{s_3}');

% figure(2);
% subplot(2,1,1);
% plot(s_last(1,:),'b','LineWidth',LineWidth);
% grid on;
% xlabel('Iteration number');
% ylabel('s_1');
% 
% subplot(2,1,2);
% plot(s_last(2,:),'b','LineWidth',LineWidth);
% grid on;
% xlabel('Iteration number');
% ylabel('s_2');

% s_1 = s(1,:,:);
% [xx,yy] = meshgrid(1:nIter,t(1):Ts:t(end));
% plot3(xx(:),yy(:),s_1(:),'b');
% grid on;
% xlabel('Iteration number');
% ylabel('Time (sec)');
% zlabel('Sliding surface');
% 
% suptitle('Iterative sliding surfaces');

%% Plant inputs in last iteration number
figure(3);
subplot(3,1,1);
plot(t(20:end),u(1,20:end,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\tau_1 (N.m)');

title('Plant inputs in last iteration number');

subplot(3,1,2);
plot(t(20:end),u(2,20:end,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\tau_2 (N.m)');

subplot(3,1,3);
plot(t(20:end),u(3,20:end,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\tau_3 (N.m)');

%% Tracking error of joints' angular velocities in different iteration numbers
figure(4);
subplot(2,1,1);
plot(t,e(1,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_1 (rad/sec)');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,e(1,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_1 (rad/sec)');
title(['Iteration #',num2str(nIter)]);

% suptitle('Tracking error of angular velocity of joint #1 in different iteration numbers');

figure(5);
subplot(2,1,1);
plot(t,e(2,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_2 (rad/sec)');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,e(2,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_2 (rad/sec)');
title(['Iteration #',num2str(nIter)]);

% suptitle('Tracking error of angular velocity of joint #2 in different iteration numbers');

figure(6);
subplot(2,1,1);
plot(t,e(3,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_3 (rad/sec)');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,e(3,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('e_3 (rad/sec)');
title(['Iteration #',num2str(nIter)]);

% suptitle('Tracking error of angular velocity of joint #3 in different iteration numbers');

%% Joints' angular velocities in different iteration numbers
figure(7);
subplot(2,1,1);
plot(t,y_d(1,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(1,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_1 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,y_d(1,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(1,:,nIter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_1 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(nIter)]);

% suptitle('Angular velocity of joint #1 in different iteration numbers');

figure(8);
subplot(2,1,1);
plot(t,y_d(2,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(2,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_2 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,y_d(2,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(2,:,nIter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_2 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(nIter)]);

% suptitle('Angular velocity of joint #2 in different iteration numbers');

figure(9);
subplot(2,1,1);
plot(t,y_d(3,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(3,:,my_iter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_3 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(my_iter)]);

subplot(2,1,2);
plot(t,y_d(3,:),'r--','LineWidth',LineWidth_d);
hold on;
plot(t,y(3,:,nIter),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('qdot_3 (rad/sec)');
legend('Desired','Real');
title(['Iteration #',num2str(nIter)]);

% suptitle('Angular velocity of joint #3 in different iteration numbers');

%% IAE of tracking errors of joints' angular velocities in all iterations
figure(10);
subplot(3,1,1);
plot(e_IAE(1,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Iteration number');
ylabel('Learning error of joint #1');

title('IAE of tracking errors in different iteration numbers');

subplot(3,1,2);
plot(e_IAE(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Iteration number');
ylabel('Learning error of joint #2');

subplot(3,1,3);
plot(e_IAE(3,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Iteration number');
ylabel('Learning error of joint #3');
