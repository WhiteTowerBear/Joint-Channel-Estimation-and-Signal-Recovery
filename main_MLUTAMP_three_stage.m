%========================================================a
%  Author: Wei Li;
%  Date:   2021.Jan.19,  SUTD;
%  Version: V1.0 
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
beta=0.2; % sparsity of transmitted signal

M=200; % The number of antennas at BS
K=500; % The number of users
T=300; % The numebr of time slots
T_p=150;% The number of pilot symbols
N=400; % The number of elements placed on RIS
K_p=200; % The number of pilot in channel H^r
N_p=0; % The number of pilot in channel H^b
damping=1;
 
var_channel=1;    % The variance of reflecting channel 

FRAME=10;
iter=30;
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
    X=sqrt(var_channel)*(randn(M,T) );
    S_loc=binornd(1,beta,M,T); 
    X=X.*S_loc;  
    X_p=X(:,1:T_p);
        
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel
    Hb=sqrt(var_channel)*(randn(N,M) );
    Hb_loc=binornd(1,beta,N,M); 
%     Hb_loc=ones(N,M);
    Hb=Hb.*Hb_loc;
 
     Hr= sqrt(var_channel)*(randn(K,N) ); 
%     [cor_U,~,~]=svd(randn(K,K));
    Hr_loc=binornd(1,beta,K,N);
%     Hr_loc=ones(K,N);
    Hr=Hr.*Hr_loc;   
%     Hr=cor_U*Hr;
%     Hr_loc=ones(K,N);
    Hr_p=Hr(1:K_p,:);
    
    % Received signals
    noise=sqrt(var_noise)*(randn(K,T));
    Rec_Y=Hr*Hb*X+noise;
    
    Y1_p=Rec_Y(1:K_p,1:T_p)*pinv(X_p);
    Y2_p=Rec_Y(1:K_p,:);
    
    Rec_Y_pilot=Rec_Y(1:K_p,:);
    TEMP=Rec_Y';
    Rec_Y_pilot_X=TEMP(1:T_p,:);
 

    %%========================================================
   %% Receiver
    % Three stage estimation: Estimate Hb with pilots in Hr and X
    % The First Stage
    [U_A1,L_A1,V_A1]=svd(Hr_p);
    A_eq_new1=L_A1*V_A1';
    rec_new1=U_A1'*Y1_p;
    prior_H1_mean=zeros(K_p,N); 
    prior_H1_var=ones(K_p,N);
    prior_X1_mean=zeros(N,M);prior_X1_var=ones(N,M);
    Y_noise1=var_noise*ones(K_p,M);
    [b_hat,v_b,a1_hat,v_a1]=Single_GAMP(prior_X1_mean,prior_X1_var,prior_H1_mean,prior_H1_var,A_eq_new1,K_p,damping,rec_new1,Y_noise1,Hb_loc,var_channel,iter,Hb,0);
     
    % The Second Stage: Estimate signal X with estimated Hb (b_hat) and
    % pilot in Hr
    [U_A2,L_A2,V_A2]=svd(Hr_p*b_hat);
     A_eq_new2=L_A2*V_A2';
     rec_new2=U_A2'*Y2_p;
     prior_H2_mean=zeros(K_p,M); % priors of Hr
    prior_H2_var=ones(K_p,M);
    % priors of X: the input X of each layer is different 
    prior_X2_mean=zeros(M,T);prior_X2_var=ones(M,T);
    Y_noise2=var_noise*ones(K_p,T);
     [x_hat,v_x,q2_hat,v_q2]=Single_GAMP(prior_X2_mean,prior_X2_var,prior_H2_mean,prior_H2_var,A_eq_new2,K_p,damping,rec_new2,Y_noise2,S_loc,var_channel,iter,X,T_p);
      
    % The Third Stage: Estimate Hr with estimated Hb (b_hat) and estimated X (x_hat)
    [U_A3,L_A3,V_A3]=svd((b_hat*x_hat)');
     A_eq_new3=L_A3*V_A3';
     rec_new3=U_A3'*Rec_Y';
     prior_H3_mean=zeros(T,N);  
     prior_H3_var=ones(T,N);
     prior_X3_mean=zeros(N,K);prior_X3_var=ones(N,K);
     Y_noise3=var_noise*ones(T,K);
     [r_trans_hat,v_r_trans,q3_hat,v_q3]=Single_GAMP(prior_X3_mean,prior_X3_var,prior_H3_mean,prior_H3_var,A_eq_new3,T,damping,rec_new3,Y_noise3,Hr_loc',var_channel,iter,Hr',0);
     r_hat=r_trans_hat';
     
     for tt=1:iter
         % The First Stage
        [U_A1,L_A1,V_A1]=svd(Hr_p);
        A_eq_new1=L_A1*V_A1';
        rec_new1=U_A1'*Rec_Y*pinv(X);
        prior_H1_mean=zeros(K,N); 
        prior_H1_var=ones(K,N);
        prior_X1_mean=zeros(N,M);prior_X1_var=ones(N,M);
        Y_noise1=var_noise*ones(K,M);
        [b_hat,v_b,a1_hat,v_a1]=Single_GAMP(prior_X1_mean,prior_X1_var,prior_H1_mean,prior_H1_var,A_eq_new1,K,damping,rec_new1,Y_noise1,Hb_loc,var_channel,iter,Hb,0);
         
        % The Second Stage: Estimate signal X with estimated Hb (b_hat) and
        % pilot in Hr
        [U_A2,L_A2,V_A2]=svd(Hr_p*b_hat);
        A_eq_new2=L_A2*V_A2';
        rec_new2=U_A2'*Y2_p;
        prior_H2_mean=zeros(K_p,M); % priors of Hr
        prior_H2_var=ones(K_p,M);
        % priors of X: the input X of each layer is different 
        prior_X2_mean=zeros(M,T);prior_X2_var=ones(M,T);
        Y_noise2=var_noise*ones(K_p,T);
        [x_hat,v_x,q2_hat,v_q2]=Single_GAMP(prior_X2_mean,prior_X2_var,prior_H2_mean,prior_H2_var,A_eq_new2,K_p,damping,rec_new2,Y_noise2,S_loc,var_channel,iter,X,T_p);

     end
      %%=================================================
%%
% Simulation result
      MSE_X=norm(X- x_hat,'fro')^2/(norm(X,'fro')^2);
      MSE_Hb=norm(Hb -b_hat,'fro')^2/(norm(Hb,'fro')^2);
      MSE_Hr=norm(Hr- r_hat,'fro')^2/(norm(Hr,'fro')^2);
%       MSE_part=norm(Hb*X -b_hat*x_hat,'fro')^2/(norm(Hb*x_hat,'fro')^2);
%       MSE_part2=norm(Hb*X -A_eq_hat_new.','fro')^2/(norm(Hb*x_hat,'fro')^2);
      MSE_all=norm(Hr*Hb*X -r_hat*b_hat*x_hat,'fro')^2/(norm(Hr*Hb*x_hat,'fro')^2);
      mc_sum_X=mc_sum_X+MSE_X;
      mc_sum_Hb=mc_sum_Hb+MSE_Hb;
      mc_sum_Hr=mc_sum_Hr+MSE_Hr;
      mc_sum_all=mc_sum_all+MSE_all;
%%========================================================%
    end
    diff_X(1,snr)=mc_sum_X/FRAME;
    diff_Hb(1,snr)=mc_sum_Hb/FRAME;
    diff_Hr(1,snr)=mc_sum_Hr/FRAME;
    diff_all(1,snr)=mc_sum_all/FRAME;

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
 
title(['UTAMP+3 steps, M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',Tp=',num2str(T_p),',Kp=',num2str(K_p)]);
xlabel('SNR (dB)');
ylabel('NMSE');
grid on
 
