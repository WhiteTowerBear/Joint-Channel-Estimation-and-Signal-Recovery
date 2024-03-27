function [x_hat,v_x,q_hat,v_q]=Single_GAMP(prior_X_mean,prior_X_var,prior_H_mean,prior_H_var,Hr,K_p,damping,Rec_Y,Y_noise,S_loc,var_channel,iter,X,T_p)
% This function is BiGAMP 
% the observation is K*N, and the estimated X is N*T, the received is K*T
 % ; part of observation Hr(1:K_p,:) is known, and the sparsity of X
    % is also known
    %The first T_p columns of X is assumed to be known
% Initialization
[K,N]=size(prior_H_mean);
[~,T]=size(prior_X_mean);

    q_hat=sqrt(var_channel)*(randn(K,N) );
    q_hat(1:K_p,:)=Hr(1:K_p,:); % The known part
    v_q=ones(K,N);
    v_q(1:K_p,:)=0;
%  
    x_hat=sqrt(var_channel)*(randn(N,T) );
    x_hat(:,1:T_p)=X(:,1:T_p);
    v_x=ones(N,T); 
    v_x(1,1:T_p)=0;
    
    q_hat_former=zeros(K,N);
 
    s_tilde_kt=zeros(K,T);
    vs_tilde_kt=ones(K,T);
 
 
    V_kt=ones(K,T);
    Z_kt=zeros(K,T);
    s_tilde_kt_former=zeros(K,T);
    vs_tilde_kt_former=ones(K,T);

    V_kt_former=ones(K,T);
 
    
    % Iterations
    diff=zeros(iter,1); 
    for tt=1:iter 
        if any(tt==[1:3])
            [Sigma_x,R_x,Sigma_b(1:K_p,:),R_b(1:K_p,:),s_tilde_kt(1:K_p,:),vs_tilde_kt(1:K_p,:)]=Single_Back_iter(Rec_Y(1:K_p,:),Y_noise(1:K_p,:),Z_kt(1:K_p,:),V_kt(1:K_p,:),q_hat(1:K_p,:),x_hat,v_q(1:K_p,:),v_x,s_tilde_kt_former(1:K_p,:),vs_tilde_kt_former(1:K_p,:),damping,tt);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
  
           [Z_kt(1:K_p,:),V_kt(1:K_p,:),V_bar_kt(1:K_p,:),q_hat(1:K_p,:),x_hat,v_q(1:K_p,:),v_x]=Single_Forw_iter(prior_X_mean,prior_X_var,prior_H_mean(1:K_p,:),prior_H_var(1:K_p,:),Sigma_x,R_x,Sigma_b(1:K_p,:),R_b(1:K_p,:),s_tilde_kt(1:K_p,:),Hr,K_p,V_kt_former(1:K_p,:),q_hat_former(1:K_p,:),damping,tt,S_loc,X,T_p);
           V_kt_former=V_kt;
           q_hat_former=q_hat;
        else
            [Sigma_x,R_x,Sigma_b,R_b,s_tilde_kt,vs_tilde_kt]=Single_Back_iter(Rec_Y,Y_noise,Z_kt,V_kt,q_hat,x_hat,v_q,v_x,s_tilde_kt_former,vs_tilde_kt_former,damping,tt);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
 
            x_dec=x_hat; % the joint estimation
        
           [Z_kt,V_kt,V_bar_kt,q_hat,x_hat,v_q,v_x]=Single_Forw_iter(prior_X_mean,prior_X_var,prior_H_mean,prior_H_var,Sigma_x,R_x,Sigma_b,R_b,s_tilde_kt,Hr,K_p,V_kt_former,q_hat_former,damping,tt,S_loc,X,T_p);
           V_kt_former=V_kt;
           q_hat_former=q_hat;
        
            diff(tt,1)=norm(x_dec-x_hat,'fro')^2/norm(x_hat,'fro')^2;
            if diff(tt,1)<1e-5
                break;
            end
        end
 
    end
%     
%     scale=X(1,:)./x_hat(1,:);
%     x_hat=x_hat.*repmat(scale,N,1);
%     q_hat=q_hat./repmat(scale.',1,N);
end