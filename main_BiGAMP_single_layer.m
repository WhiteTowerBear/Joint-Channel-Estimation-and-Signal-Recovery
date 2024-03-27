%========================================================a
%  Author: Wei Li;
%  Date:   2021.Jan.19,  SUTD;
%  Version: V1.0 
%  Note: This is the code for Single layer BiGAMP, including
%        BiGAMP and UTAMP
%%========================================================
clear all
clc

fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
global S_loc
beta=0.2; % sparsity of transmitted signal

M=10; % The number of antennas at BS
K=32; % The number of users
T=48; % The numebr of time slots
T_p=0;% The number of pilot symbols
N=40; % The number of elements placed on RIS
K_p=20; % The number of pilot in channel H^r
N_p=0; % The number of pilot in channel H^b
damping=0.7;
 
var_channel=1;    % The variance of reflecting channel 

FRAME=10;
iter=50;
diff_temp=[];
snr=1;


for SNR=0:5:30  
    SNR
    mc_sum_Hb=0;
    mc_sum_Hr=0;
    mc_sum_X=0;
    mc_sum_all=0; % the joint estimation
    
    for frame=1:1:FRAME
    % Generate transmitted signals X
    
    tranx=sqrt(var_channel)*(randn(N,T) );
    S_loc=binornd(1,beta,N,T); 
    X=tranx.*S_loc;

%     tranx=sqrt(var_channel)*(randn(M,T) );
%     S_loc=binornd(1,beta,M,T); 
%     tranx=tranx.*S_loc;  
%  
%     Hb=sqrt(var_channel)*(randn(N,M) );
%     Hb_loc=binornd(1,beta,N,M); 
%     Hb=Hb.*Hb_loc;
%     
%     X=Hb*tranx;
%     HbX_loc=Hb_loc*S_loc;
%     HbX_loc=1-(HbX_loc==0);
  
        
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel
    Hr= sqrt(var_channel)*(randn(K,N) ); 
    S_loc_HR=ones(K,N);
%     S_loc_HR=binornd(1,beta,K,N); 
    Hr=Hr.*S_loc_HR;
    % Received signals
    noise=sqrt(var_noise)*(randn(K,T));
    Rec_Y=Hr*X+noise;
    
    
    
 
    %%========================================================
   %% Receiver
    % Backward message passing
    prior_H_mean=zeros(K,N); % priors of Hr
    prior_H_var=ones(K,N);
    % priors of X: the input X of each layer is different 
    prior_X_mean=zeros(N,T);prior_X_var=ones(N,T);
    Y_noise=var_noise*ones(K,T);
    % the observation is K*N, and the estimated X is N*T, the received is
    % K*T ; part of observation Hr(1:K_p,:) is known, and the sparsity of X
    % is also known
    
    method='BiGAMP';
    switch(method)
        case 'BiGAMP'
            [x_hat,v_x,q_hat,v_q]=Single_GAMP(prior_X_mean,prior_X_var,prior_H_mean,prior_H_var,Hr,K_p,damping,Rec_Y,Y_noise,S_loc,var_channel,iter,X,T_p);
        case 'UTAMP'
        % Unitary Transformation AMP (UTAMP)
             [U_A,L_A,V_A]=svd(Hr);
             A_eq_new=L_A*V_A';
             rec_new=U_A'*Rec_Y;
             [x_hat,v_x,q_hat,v_q]=Single_GAMP(prior_X_mean,prior_X_var,prior_H_mean,prior_H_var,A_eq_new,K,damping,rec_new,Y_noise,S_loc,var_channel,iter,X,T_p);
             q_hat=U_A*q_hat;
    end
 
    %%=================================================
%%
% Simulation result
      MSE_X=norm(X- x_hat,'fro')^2/(norm(X,'fro')^2);
      MSE_Hr=norm(Hr- q_hat,'fro')^2/(norm(Hr,'fro')^2);
      MSE_all=norm(Hr*X -q_hat*x_hat,'fro')^2/(norm(Hr*X,'fro')^2);
      mc_sum_X=mc_sum_X+MSE_X;
 
      mc_sum_Hr=mc_sum_Hr+MSE_Hr;
      mc_sum_all=mc_sum_all+MSE_all;
%%========================================================%
    end
    diff_X(1,snr)=mc_sum_X/FRAME;
 
    diff_Hr(1,snr)=mc_sum_Hr/FRAME;
    diff_all(1,snr)=mc_sum_all/FRAME;

    snr=snr+1;
end
fprintf(fid,'\n \n M=%d  K=%d  N=%d  T=%d  %d Monte Carlo simulations \n',M,K,N,T,FRAME);
fprintf(fid,' %s\r\n ',datestr(now,31));

 
fprintf(fid,'NMSE of Hr is \n');
fprintf(fid, ' %d  ', diff_Hr);
fprintf(fid,'\n NMSE of X is \n');
fprintf(fid, ' %d  ', diff_X);
fclose(fid);
figure

semilogy(0:5:30,diff_X,'-ro','linewidth',2)
hold on
 
semilogy(0:5:30,diff_Hr,'-bd','linewidth',2)
hold on
semilogy(0:5:30,diff_all,'-k+','linewidth',2)
hold on
legend({'Proposed MP (X)', 'Proposed MP (H^r)','Proposed MP (all)'});
 
% hold on; 
% h = zeros(4, 1); 
% h(1) = plot(NaN,NaN,'d:k','linewidth',2); 
% h(2) = plot(NaN,NaN,'s-k','linewidth',2); 
% h(3) = plot(NaN,NaN,'o:k','linewidth',2); 
% h(4) = plot(NaN,NaN,'x-k','linewidth',2);
% legend(h, 'MRT (IP)','MRT (P)','ZF (IP)','ZF (P)');
switch(method)
    case 'BiGAMP'
        title(['BiGAMP,', 'M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',\beta=',num2str(damping),',Tp=',num2str(T_p),',Kp=',num2str(K_p),',Np=',num2str(N_p)]);
    case 'UTAMP'
        title(['UTAMP,', 'M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',\beta=',num2str(damping),',Tp=',num2str(T_p),',Kp=',num2str(K_p),',Np=',num2str(N_p)]);
end
xlabel('SNR (dB)');
ylabel('NMSE');
grid on
 
