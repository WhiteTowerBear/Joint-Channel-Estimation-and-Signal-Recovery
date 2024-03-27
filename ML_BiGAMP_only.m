function [q_hat,b_hat,x_hat]=ML_BiGAMP_only(Rec_Y,Y_noise,Hr,Hb,X,var_channel,Hr_loc,Hb_loc,S_loc,HbX_loc,K_p,T_p,N_p,iter,damping,damping_r,Rec_Y_pilot)
% This function is the part of two-layers BiGAMP algorithm
[K,T]=size(Rec_Y);
[N,M]=size(Hb);

    prior_H=[0,var_channel]; % priors of H^b
    prior_Q=[0,var_channel]; % priors of Q
    % priors of X: the input X of each layer is different 
    prior_X_mean=zeros(M,T);prior_X_var=ones(M,T);
     % Initialization
    
    q_hat=zeros(K,N);
    q_hat=q_hat.*Hr_loc;
    q_hat(1:K_p,:)=Hr(1:K_p,:); % The known part
    v_q=ones(K,N);
    v_q=v_q.*Hr_loc;
    v_q(1:K_p,:)=0;
     
    b_hat=sqrt(var_channel)*(randn(N,M) );  
    b_hat(1:N_p,:)=Hb(1:N_p,:);
    b_hat=b_hat.*Hb_loc;
    v_b=ones(N,M);
    v_b(1:N_p,:)=v_b(1:N_p,:);
    v_b=v_b.*Hb_loc;
 
    
    x_hat=sqrt(var_channel)*(randn(M,T) );
    x_hat=x_hat.*S_loc; 
    x_hat(:,1:T_p)=X(:,1:T_p); % The pilot data part
    v_x=ones(M,T);
    v_x=v_x.*S_loc; 
    v_x(1,1:T_p)=0;
    
    x_hat_former=zeros(M,T);b_hat_former=zeros(N,M);
    u_hat_former=zeros(N,T);q_hat_former=zeros(K,N);
    
    u_hat=b_hat*x_hat;
    v_u=ones(N,T);
    Z_nt=zeros(N,T);V_nt=ones(N,T); 
    Z_kt=zeros(K,T);V_kt=ones(K,T);
    Z_nt_pilot=zeros(N_p,T_p);
    V_nt_pilot=zeros(N_p,T_p);
    if (T_p~=0 && N_p~=0)
        Z_nt_pilot=Hb(1:N_p,:)*X(:,1:T_p);
        V_nt_pilot=zeros(N_p,T_p);
        Z_nt(1:N_p,1:T_p)=Z_nt_pilot;
        V_nt(1:N_p,1:T_p)=0;
    end
    
    s_tilde_kt=zeros(K,T);
    vs_tilde_kt=ones(K,T);
    s_tilde_nt=zeros(N,T);
    vs_tilde_nt=ones(N,T);
    V_nt=ones(N,T);
    V_kt=ones(K,T);
    Z_nt=zeros(N,T);
    Z_kt=zeros(K,T);
    s_tilde_kt_former=zeros(K,T);
    vs_tilde_kt_former=ones(K,T);
    s_tilde_nt_former=zeros(N,T);
    vs_tilde_nt_former=ones(N,T);
    V_nt_former=ones(N,T);
    V_kt_former=ones(K,T);
    Z_nt_former=zeros(N,T);
    Z_kt_former=zeros(K,T);
    V_bar_nt_former=zeros(N,T);
    V_bar_kt_former=zeros(K,T);
    V_bar_kt=zeros(K,T);
    
    % Iterations
    x_dec=zeros(M,T);
  
    
        diff=ones(iter,1);
        for tt=1:iter
        % Hr part
        if any(tt==[1:2])
            [Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt(1:K_p,1:T_p),vs_tilde_kt(1:K_p,1:T_p)]=Back_iter(2,Rec_Y_pilot(:,1:T_p),Y_noise(1:K_p,1:T_p),Z_kt(1:K_p,1:T_p),V_kt(1:K_p,1:T_p),q_hat(1:K_p,:),u_hat(:,1:T_p),v_q(1:K_p,:),v_u(:,1:T_p),s_tilde_kt_former(1:K_p,1:T_p),vs_tilde_kt_former(1:K_p,1:T_p),T_p,N_p,damping,tt,HbX_loc(:,1:T_p),Hb,X);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
            
            [Sigma_x(:,1:T_p),R_x(:,1:T_p),Sigma_b,R_b,s_tilde_nt(:,1:T_p),vs_tilde_nt(:,1:T_p)]=Back_iter(1,R_u(:,1:T_p),Sigma_u(:,1:T_p),Z_nt(:,1:T_p),V_nt(:,1:T_p),b_hat,x_hat(:,1:T_p),v_b,v_x(:,1:T_p),s_tilde_nt_former(:,1:T_p),vs_tilde_nt_former(:,1:T_p),T_p,N_p,damping,tt,HbX_loc(:,1:T_p),Hb,X);
            s_tilde_nt_former=s_tilde_nt;
            vs_tilde_nt_former=vs_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
            
        
           [Z_nt(:,1:T_p),V_nt(:,1:T_p),V_bar_nt(:,1:T_p),b_hat,x_hat(:,1:T_p),v_b,v_x(:,1:T_p)]=Forw_iter(1,prior_X_mean(:,1:T_p),prior_X_var(:,1:T_p),prior_H,Sigma_x(:,1:T_p),R_x(:,1:T_p),Sigma_b,R_b,s_tilde_nt(:,1:T_p),X,T_p,Hr(1:K_p,:),K_p,V_nt_former(:,1:T_p),b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc(:,1:T_p),Hb_loc,HbX_loc(:,1:T_p),Hr_loc(1:K_p,:));
           V_nt_former=V_nt;
           b_hat_former=b_hat;
           [Z_kt(1:K_p,1:T_p),V_kt(1:K_p,1:T_p),V_bar_kt(1:K_p,1:T_p),q_hat(1:K_p,:),u_hat(:,1:T_p),v_q(1:K_p,:),v_u(:,1:T_p)]=Forw_iter(2,Z_nt(:,1:T_p),V_nt(:,1:T_p),prior_Q,Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,1:T_p),X,T_p,Hr(1:K_p,:),K_p,V_kt_former(1:K_p,1:T_p),q_hat_former(1:K_p,:),N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc(:,1:T_p),Hb_loc,HbX_loc(:,1:T_p),Hr_loc(1:K_p,:)); 
           V_kt_former=V_kt; 
           q_hat_former=q_hat;
             
        else
            [Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,:),vs_tilde_kt(1:K_p,:)]=Back_iter(2,Rec_Y_pilot,Y_noise(1:K_p,:),Z_kt(1:K_p,:),V_kt(1:K_p,:),q_hat(1:K_p,:),u_hat,v_q(1:K_p,:),v_u,s_tilde_kt_former(1:K_p,:),vs_tilde_kt_former(1:K_p,:),T_p,N_p,damping,tt,HbX_loc,Hb,X);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
            [Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,vs_tilde_nt]=Back_iter(1,R_u,Sigma_u,Z_nt,V_nt,b_hat,x_hat,v_b,v_x,s_tilde_nt_former,vs_tilde_nt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X);
            s_tilde_nt_former=s_tilde_nt;
            vs_tilde_nt_former=vs_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
            
        
           [Z_nt,V_nt,V_bar_nt,b_hat,x_hat,v_b,v_x]=Forw_iter(1,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,X,T_p,Hr(1:K_p,:),K_p,V_nt_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc(1:K_p,:));
           V_nt_former=V_nt;
           b_hat_former=b_hat;
           [Z_kt(1:K_p,:),V_kt(1:K_p,:),V_bar_kt(1:K_p,:),q_hat(1:K_p,:),u_hat,v_q(1:K_p,:),v_u]=Forw_iter(2,Z_nt,V_nt,prior_Q,Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,:),X,T_p,Hr(1:K_p,:),K_p,V_kt_former(1:K_p,:),q_hat_former(1:K_p,:),N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc(1:K_p,:)); 
           V_kt_former=V_kt; 
           q_hat_former=q_hat;
           
 
            diff(tt,1)=norm(x_dec-b_hat*x_hat,'fro')^2/norm(b_hat*x_hat,'fro')^2;
           if diff(tt,1)<1e-5
               break;
           end           
        end
        end %  outer loop stops (pilot parts) 
        

%         
        for tt=1:10
             [Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt,vs_tilde_kt]=Back_iter(2,Rec_Y,Y_noise,Z_kt,V_kt,q_hat,u_hat,v_q,v_u,s_tilde_kt_former,vs_tilde_kt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
            [Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,vs_tilde_nt]=Back_iter(1,R_u,Sigma_u,Z_nt,V_nt,b_hat,x_hat,v_b,v_x,s_tilde_nt_former,vs_tilde_nt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X);
            s_tilde_nt_former=s_tilde_nt;
            vs_tilde_nt_former=vs_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
        
           [Z_nt,V_nt,V_bar_nt,b_hat,x_hat,v_b,v_x]=Forw_iter(1,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,X,T_p,q_hat,K,V_nt_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,0.7,tt,S_loc,Hb_loc,HbX_loc,Hr_loc);
           V_nt_former=V_nt;
           b_hat_former=b_hat;
           [Z_kt,V_kt,V_bar_kt,q_hat,u_hat,v_q,v_u]=Forw_iter(2,Z_nt,V_nt,prior_Q,Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt,X,T_p,q_hat,K,V_kt_former,q_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,0.7,tt,S_loc,Hb_loc,HbX_loc,Hr_loc); 
           V_kt_former=V_kt; 
           q_hat_former=q_hat;
           
           diff1(tt,1)=norm(x_dec-b_hat*x_hat,'fro')^2/norm(b_hat*x_hat,'fro')^2;
        end
         
end