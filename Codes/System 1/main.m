clc;
clear;
close all;

tic;
%% Time sequence
T0 = 0;
Ts = 0.001;
Tf = 18;

t = T0:Ts:Tf;

nSteps = length(t);

%% Controller parameters
% MFAILITSMC
eta = 1;
mu = 1;

epsilon = 1e-5;

rho = 3/5;

c_1 = 10*diag([1,1]);
c_2 = 5*diag([0.1,0.1]);

c_s = 0.0001*diag([0.002,0.002]);

m_1 = 0*diag([1,1]);
m_2 = 0.1*diag([0.1,1]);
m_3 = 0*[1 1]';
m = [m_1 m_2 m_3];

nIter = 100;

%% Initialization
nInputs = 2;
nOutputs = 2;

u = zeros(nInputs,nSteps,nIter);
u(:,:,1) = zeros(nInputs,nSteps);
delta_u_eq = zeros(nInputs,nSteps,nIter);
delta_u_sw = zeros(nInputs,nSteps,nIter);
delta_u = zeros(nInputs,nSteps,nIter);

y = zeros(nOutputs,nSteps+1,nIter);
y(:,:,1) = zeros(nOutputs,nSteps+1);
delta_y = zeros(nOutputs,nSteps+1,nIter);

delta_omega = zeros(nOutputs+nInputs+1,nSteps+1,nIter);

y_d = zeros(nOutputs,nSteps+1);
y_d(:,1) = [ (9*pi/40)*cos((9*pi/40)*t(1))
             (7*pi/16)*cos((7*pi/32)*t(1)) ];
y_d(:,1) = [ sat_func_yd(y_d(1,1),-0.5,0.5)
             sat_func_yd(y_d(2,1),-0.5,0.5) ];

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
e_rms(:,1) = 3*[1 1]';
e_mean = zeros(nOutputs,nIter);
e_mean(:,1) = [0.006 0.002]';
e_IAE = zeros(nOutputs,nIter);
e_IAE(:,1) = [1.3160 2.9161]';

fprintf('Simulating MFAILITSMC on 2-DOF exoskeleton...\n');
fprintf('Number of iterations is %d\n',nIter);
fprintf('-----------------------------------------------\n');

%% Iteration loop
for i = 2:nIter
    
    fprintf('Iteration #%d is evaluated...\n',i);
    
    q = [0 0]';
    q_dot = 0.01*rand(2,1);
%     q_dot = [0 0]';
    y(:,1,i) = q_dot;

    for k = 1:nSteps
        
        %% PJM estimation
        if(i == 2)
            
            PHI_hat_1{k,1} = 1*diag([1,1]);
            PHI_hat_2{k,1} = 0.5*diag([.005,.001]);
            PHI_hat_3{k,1} = 1*[0 0]';
            
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
%         y_d(:,k+1) = [ (3*pi/25)*cos((pi/25)*(k+1)*Ts)
%                         (pi/10)*cos((pi/50)*(k+1)*Ts)  ];
                    
%         y_d(:,k+1) = [ 1
%                        1 ];

%         y_d(:,k+1) = [ 1.5*cos((2*pi/4)*t(k))
%                        1.2*cos((2*pi/7)*t(k)) ];

        y_d(:,k+1) = [ (9*pi/40)*cos((9*pi/40)*t(k))
                       (7*pi/16)*cos((7*pi/32)*t(k)) ];

        y_d(:,k+1) = [ sat_func_yd(y_d(1,k+1),-0.5,0.5)
                       sat_func_yd(y_d(2,k+1),-0.5,0.5) ];
    
        %% MFAILITSMC
        delta_u_eq(:,k,i) = ((PHI_hat_2{k,i}+m_2)^-1)*(y_d(:,k+1)-y(:,k+1,i-1)- ...
                            (PHI_hat_1{k,i}+m_1)*delta_y(:,k,i)-PHI_hat_3{k,i}+...
                            (c_1^-1)*c_2*int_e(:,k+1,i-1)-(c_1^-1)*s(:,k+1,i-1));
        
        delta_u_sw(:,k,i) = c_s*((PHI_hat_2{k,i}+m_2)^-1)*sign(s(:,k,i));
%         delta_u_sw(:,k,i) = c_s*((PHI_hat_2{k,i}+m_2)^-1)*tanh(s(:,k,i))/0.1;

        delta_u(:,k,i) = delta_u_eq(:,k,i) + delta_u_sw(:,k,i);
        u(:,k,i) = u(:,k,i-1) + delta_u(:,k,i);

        %% Plant simulation
%         dist(:,k) = [   1*exp(i)*(t(k)>=5 & t(k)<=13) 
%                       0.5*exp(i)*(t(k)>=7 & t(k)<=15) ];

%         dist(:,k) = [ 3*sin(2*t(k))+0.1*cos(i*pi/3) 
%                       2*cos(1*t(k))+0.2*sin(i*pi/5) ];

%         dist(:,k) = [ 5*sin(2*t(k)) 
%                       2*cos(1*t(k)) ];
                  
        dist(:,k) = [ 5*sin(2*t(k)) + 1*0.002*cos(2*pi/i)
                      2*cos(1*t(k)) + 1*0.002*sin(3*pi/i) ];
                  
%         w(2,k,i) = (1+0.1*cos(i*pi/15))*exp(-0.05*k) + (0.2+0.1*sin(pi/i))*sin(2*k);
%         w(3,k,i) = 0.5 + 0.1*sin(i*pi/15) + 0.1*cos(k/5);
%         w(4,k,i) = 0.5 + 0.2*cos(2*pi/i)*((-1)^(-5*i*k));
                  
        q_ddot = twoDOF_exoskeleton_dynamics_imp(q,q_dot,u(:,k,i),dist(:,k));
        q_dot = q_dot + Ts*q_ddot;
        q = q + Ts*q_dot;
        
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

fprintf('\nSuccessfully done!\n');
fprintf('-----------------------------------------------\n');

toc;
%% Plot results
LineWidth = 1.5;
LineWidth_d = 2;

my_iter = 10;

%% Difference of iterative sliding surfaces in last iteration number
figure(1);
subplot(2,1,1);
plot(t,delta_s(1,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\Delta{s_1}');

title('Difference of sliding surfaces in last iteration number');

subplot(2,1,2);
plot(t,delta_s(2,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\Delta{s_2}');

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
subplot(2,1,1);
plot(t,u(1,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\tau_1 (N.m)');

title('Plant inputs in last iteration number');

subplot(2,1,2);
plot(t,u(2,:,end),'b','LineWidth',LineWidth);
grid on;
xlabel('Time (sec)');
ylabel('\tau_2 (N.m)');

%% Tracking errors of joints' angular velocities in different iteration numbers
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

%% Joints' angular velocities in different iteration numbers
figure(6);
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

figure(7);
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

%% IAE of tracking errors of joints' angular velocities in all iterations
figure(8);
subplot(2,1,1);
plot(e_IAE(1,1:50),'b','LineWidth',LineWidth);
grid on;
xlabel('Iteration number');
ylabel('Learning error of joint #1');

title('IAE of tracking errors in different iteration numbers');

subplot(2,1,2);
plot(e_IAE(2,:),'b','LineWidth',LineWidth);
grid on;
xlabel('Iteration number');
ylabel('Learning error of joint #2');
