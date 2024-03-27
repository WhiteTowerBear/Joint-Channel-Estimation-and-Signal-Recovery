%========================================================a
%  Author: Wei Li;
%  Date:   2021.March.24,  SUTD;
%  Version: V1.0 
%  Note: This is the code for benchmark: BiGAMP+LS
% BiGAMP estimates Hr and Hb, LS estimates X
%%========================================================
clear all
clc

fid=fopen('result.txt','a+');
%%========================================================
%  Set the parameters
global S_loc
beta=0.2; % sparsity of transmitted signal

M=100; % The number of antennas at BS
K=500; % The number of users
T=200; % The numebr of time slots
T_p=100;% The number of pilot symbols
N=200; % The number of elements placed on RIS
K_p=150; % The number of pilot in channel H^r
N_p=0; % The number of pilot in channel H^b
M_p=0;
damping=0.7;
 
var_channel=1;    % The variance of reflecting channel 

FRAME=500;
iter=50;
diff_temp=[];
snr=1;


for SNR=0:5:30  
    SNR
    mc_sum_Hb=0;
    mc_sum_Hr=0;
    mc_sum_Hb=0;
    mc_sum_X=0; % the joint estimation
    
    for frame=1:1:FRAME
     % Generate transmitted signals X
    X=sqrt(var_channel)*(randn(M,T) );
    S_loc=binornd(1,beta,M,T); 
    X=X.*S_loc;  
        
    var_noise=10^(-0.1*SNR);
    %%========================================================
    %  Pass channel

    Hb=sqrt(var_channel)*(randn(N,M) );
    Hb_loc=binornd(1,beta,N,M); 
    Hb=Hb.*Hb_loc;
    
    HbX_loc=Hb_loc*S_loc;
    HbX_loc=1-(HbX_loc==0);
     
    Hr= sqrt(var_channel)*(randn(K,N) ); 
    Hr_loc=binornd(1,beta,K,N);
    Hr=Hr.*Hr_loc;  
 
    
    Hbr_loc=Hb_loc.'*Hr_loc';
    Hbr_loc=1-(Hbr_loc==0);
    
    % Received signals
    noise=sqrt(var_noise)*(randn(K,T));
%     noise=0;
    H_eq=Hr*Hb;
    Rec_Y=H_eq*X+noise;
     
 
    %%=================================================
    % Two known factors
%     if K<N
%         Hr_inv=Hr'*inv(Hr*Hr');
%     else
%         Hr_inv=inv(Hr'*Hr)*Hr';
%     end
    Hr_inv=pinv(Hr);
%     if M<T
%         X_inv=X'*inv(X*X');
%     else
%         X_inv=inv(X'*X)*X';
%     end
    X_inv=pinv(X);
    r_mat=Hb*X;
%     if N<T
%         HbX_inv=r_mat'*inv(r_mat*r_mat');
%     else
%         HbX_inv=inv(r_mat'*r_mat)*r_mat';
%     end
    HbX_inv=pinv(r_mat);
 
%     if K<M
%         Heq_inv=H_eq'*inv(H_eq*H_eq');
%     else
%         Heq_inv=inv(H_eq'*H_eq)*H_eq';
%     end
    Heq_inv=pinv(H_eq);
    
    b_hat=Hr_inv*Rec_Y*X_inv;v_b=0;
    r_hat=Rec_Y*HbX_inv;v_r=0;
    x_hat=Heq_inv*Rec_Y;v_x=0;
    x_hat(:,1:T_p)=X(:,1:T_p);
    r_hat(1:K_p,:)=Hr(1:K_p,:);
    x_hat=x_hat.*S_loc;
    b_hat=b_hat.*Hb_loc;
    r_hat=r_hat.*Hr_loc;
%%
% Simulation result
      MSE_Hb=norm(Hb- b_hat,'fro')^2/(norm(Hb,'fro')^2);
      MSE_Hr=norm(Hr- r_hat,'fro')^2/(norm(Hr,'fro')^2);
      MSE_X=norm(X -x_hat,'fro')^2/(norm(X,'fro')^2);
      mc_sum_Hb=mc_sum_Hb+MSE_Hb;
 
      mc_sum_Hr=mc_sum_Hr+MSE_Hr;
      mc_sum_X=mc_sum_X+MSE_X;
%%========================================================%
    end
    diff_Hb(1,snr)=mc_sum_Hb/FRAME;
 
    diff_Hr(1,snr)=mc_sum_Hr/FRAME;
    diff_X(1,snr)=mc_sum_X/FRAME;

    snr=snr+1;
end
fprintf(fid,'\n \n M=%d  K=%d  N=%d  T=%d  %d Monte Carlo simulations \n',M,K,N,T,FRAME);
fprintf(fid,' %s\r\n ',datestr(now,31));

 
fprintf(fid,'NMSE of Hr is \n');
fprintf(fid, ' %d  ', diff_Hr);
fprintf(fid,'\n NMSE of X is \n');
fprintf(fid, ' %d  ', diff_Hb);
fclose(fid);
figure


semilogy(0:5:30,diff_X,'-ro','linewidth',2)
hold on
semilogy(0:5:30,diff_Hb,'-ms','linewidth',2)
hold on
semilogy(0:5:30,diff_Hr,'-bd','linewidth',2)
hold on

legend({'$\mathbf{X}$ (LS)','$\mathbf{H}^b$ (LS)','$\mathbf{H}^r$ (LS)'});
 
% hold on; 
% h = zeros(4, 1); 
% h(1) = plot(NaN,NaN,'d:k','linewidth',2); 
% h(2) = plot(NaN,NaN,'s-k','linewidth',2); 
% h(3) = plot(NaN,NaN,'o:k','linewidth',2); 
% h(4) = plot(NaN,NaN,'x-k','linewidth',2);
% legend(h, 'MRT (IP)','MRT (P)','ZF (IP)','ZF (P)');

title(['LS,', 'M=',num2str(M),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',\beta=',num2str(damping),',Tp=',num2str(T_p),',Kp=',num2str(K_p),',Np=',num2str(N_p)]);

xlabel('SNR (dB)');
ylabel('NMSE');
grid on
 
