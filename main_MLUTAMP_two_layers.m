%========================================================a
%  Author: Wei Li;
%  Date:   2021.Feb.2,  SUTD;
%  Version: V1.1 
%  Note: This is the code for MP algorithm in RIS-assisted 
%     system that estimates transmitted signal and two channels
%     without PARAFAC structure.
% Drawbacks: diverges easily with improper damping factors and sparsity
%%========================================================
clear all
clc

fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
global S_loc
global Hb_loc
global Hr_loc
iter_flag=0;

beta=0.2; % sparsity of transmitted signal

M=100; % The number of antennas at BS
K=500; % The number of users
T=200; % The numebr of time slots
T_p=100;% The number of pilot symbols
N=200; % The number of elements placed on RIS
K_p=150; % The number of pilot in channel H^r
N_p=0; % The number of pilot in channel H^b
damping=0.7;
 damping_r=1; % The damping factor in inner iteration
var_channel=1;    % The variance of reflecting channel 

FRAME=100;
iter=30;
diff_temp=[];
iter_count=[];
snr=1;


for SNR=0:5:30  
    SNR
    mc_sum_Hb=0;
    mc_sum_Hr=0;
    mc_sum_X=0;
    mc_sum_all=0; % the joint estimation
    iter_sum=0;
    
    for frame=1:1:FRAME
    % Generate transmitted signals X
%     X=sqrt(var_channel)*(randn(M,T) );
    X=randi([0,1],M,T);
    X(find(X<1))=-1;
    X=X/sqrt(2);
%     S_loc=binornd(1,beta,M,T); 
    S_loc=ones(M,T);
    X=X.*S_loc;  
    
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel

    Hb=sqrt(var_channel)*(randn(N,M) );
%     Hb_loc=binornd(1,beta,N,M); 
%     Hb=Hb.*Hb_loc;
    [Hb_loc,Hb]=correlated_matrix(Hb);
     
    HbX_loc=Hb_loc*S_loc;
    HbX_loc=1-(HbX_loc==0);
 
    Hr= sqrt(var_channel)*(randn(K,N) );
    Hr_loc=binornd(1,beta,K,N); 
    Hr=Hr.*Hr_loc;
%     [Hr_loc,Hr]=correlated_matrix(Hr);
%     Hr_loc=ones(K,N);
   
    Hbr_loc=Hb_loc.'*Hr_loc';
    Hbr_loc=1-(Hbr_loc==0);
    
    % Received signals
    noise=sqrt(var_noise)*(randn(K,T));
    Rec_Y=Hr*Hb*X+noise;
    Rec_Y_pilot=Rec_Y(1:K_p,:);
    TEMP=Rec_Y';
    Rec_Y_pilot_X=TEMP(1:T_p,:);
    Y_noise=var_noise*ones(K,T);

    %%========================================================
   %% Receiver
    % Backward message passing
    [q_hat,b_hat,x_hat,iter_flag]=ML_UTBiGAMP_EP(Rec_Y,Y_noise,Hr,Hb,X,var_channel,Hr_loc,Hb_loc,S_loc,HbX_loc,K_p,T_p,N_p,iter,damping,damping_r,Rec_Y_pilot);
%     [~,~,q_hat_new]=ML_BiGAMP_two(Rec_Y.',Y_noise.',X',Hb',Hr',var_channel,S_loc',Hb_loc',Hr_loc',Hbr_loc,T_p,K_p,M_p,iter,damping,damping_r,Rec_Y_pilot_X);
%     x_hat=x_hat_new';
%     b_hat=b_hat_new';
%     q_hat=q_hat_new';
    %%=================================================
%%
% Simulation result
      MSE_X=norm(X- x_hat,'fro')^2/(norm(X,'fro')^2);
      MSE_Hb=norm(Hb -b_hat,'fro')^2/(norm(Hb,'fro')^2);
      MSE_Hr=norm(Hr- q_hat,'fro')^2/(norm(Hr,'fro')^2);
%       MSE_part=norm(Hb*X -b_hat*x_hat,'fro')^2/(norm(Hb*x_hat,'fro')^2);
%       MSE_part2=norm(Hb*X -A_eq_hat_new.','fro')^2/(norm(Hb*x_hat,'fro')^2);
      MSE_all=norm(Hr*Hb*X -q_hat*b_hat*x_hat,'fro')^2/(norm(Hr*Hb*x_hat,'fro')^2);
      mc_sum_X=mc_sum_X+MSE_X;
      mc_sum_Hb=mc_sum_Hb+MSE_Hb;
      mc_sum_Hr=mc_sum_Hr+MSE_Hr;
      mc_sum_all=mc_sum_all+MSE_all;
      iter_sum=iter_sum+iter_flag;
%%========================================================%
    end
    diff_X(1,snr)=mc_sum_X/FRAME;
    diff_Hb(1,snr)=mc_sum_Hb/FRAME;
    diff_Hr(1,snr)=mc_sum_Hr/FRAME;
    diff_all(1,snr)=mc_sum_all/FRAME;
    iter_count(1,snr)=iter_sum/FRAME
    
    snr=snr+1;
end
fprintf(fid,'\n \n M=%d  K=%d  N=%d  T=%d  %d Monte Carlo simulations \n',M,K,N,T,FRAME);
fprintf(fid,' %s\r\n ',datestr(now,31));

fprintf(fid,'\n NMSE of Hb is \n');
fprintf(fid, ' %d  ', diff_Hb);
fprintf(fid,'NMSE of Hr is \n');
fprintf(fid, ' %d  ', diff_Hr);
fprintf(fid,'\n NMSE of X is \n');
fprintf(fid, ' %d  ', diff_X);
fclose(fid);
figure

semilogy(0:5:30,diff_X,'-ro','linewidth',2)
hold on
semilogy(0:5:30,diff_Hb,'-ms','linewidth',2)
hold on
semilogy(0:5:30,diff_Hr,'-bd','linewidth',2)
hold on
semilogy(0:5:30,diff_all,'-k+','linewidth',2)
hold on
legend({'Proposed MP (X)','Proposed MP (H^b)','Proposed MP (H^r)','Proposed MP (all)'});
 
% hold on; 
% h = zeros(4, 1); 
% h(1) = plot(NaN,NaN,'d:k','linewidth',2); 
% h(2) = plot(NaN,NaN,'s-k','linewidth',2); 
% h(3) = plot(NaN,NaN,'o:k','linewidth',2); 
% h(4) = plot(NaN,NaN,'x-k','linewidth',2);
% legend(h, 'MRT (IP)','MRT (P)','ZF (IP)','ZF (P)');
 
title(['UTAMP, M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',\beta=',num2str(damping_r),',Tp=',num2str(T_p),',Kp=',num2str(K_p),',Np=',num2str(N_p)]);
xlabel('SNR (dB)');
ylabel('NMSE');
grid on
 
