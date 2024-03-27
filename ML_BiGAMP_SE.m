function [q_hat,b_hat,x_hat,se_x,se_b,se_r]=ML_BiGAMP_SE(Rec_Y,Y_noise,Hr,Hb,X,var_channel,Hr_loc,Hb_loc,S_loc,HbX_loc,K_p,T_p,N_p,iter,damping,damping_r,Rec_Y_pilot)
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
    v_tilde_2=zeros(K,T);
    s_tilde_nt=zeros(N,T);
    vs_tilde_nt=ones(N,T);
    v_tilde_1=zeros(N,T);
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
    Vnt_scalar=1;
    Vkt_scalar_former=1;
    Vnt_scalar_former=1;
    vsnt_tilde_scalar_former=1;
    vskt_tilde_scalar_former=1;
    Vnt_scalar=1;
    Vkt_scalar=1;
    
  % State evolution
        var_noise=Y_noise(1,1);
 
    % Iterations
    x_dec=zeros(M,T);
  V1_se=1;
  V2_se=1;
  vars_z1=1;
  vars_r=1;
  vars_x=1;
  vars_b=1;  
  V1_se_former=0;
  V2_se_former=0;
  chi_x=0.2;
  chi_b=0.2;
  chi_r=0.2;
  vs2_se_former=0;
  vs1_se_former=0;
  
        diff=ones(iter,1);
        for tt=1:iter
        % Hr part
        if any(tt==[1:2])
            [Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt(1:K_p,1:T_p),vs_tilde_kt(1:K_p,1:T_p),v_tilde_2(1:K_p,1:T_p),vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(2,Rec_Y_pilot(:,1:T_p),Y_noise(1:K_p,1:T_p),Z_kt(1:K_p,1:T_p),V_kt(1:K_p,1:T_p),q_hat(1:K_p,:),u_hat(:,1:T_p),v_q(1:K_p,:),v_u(:,1:T_p),s_tilde_kt_former(1:K_p,1:T_p),vs_tilde_kt_former(1:K_p,1:T_p),T_p,N_p,damping,tt,HbX_loc(:,1:T_p),Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
 
 
            [Sigma_x(:,1:T_p),R_x(:,1:T_p),Sigma_b,R_b,s_tilde_nt(:,1:T_p),vs_tilde_nt(:,1:T_p),v_tilde_1(:,1:T_p),vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(1,R_u(:,1:T_p),Sigma_u(:,1:T_p),Z_nt(:,1:T_p),V_nt(:,1:T_p),b_hat,x_hat(:,1:T_p),v_b,v_x(:,1:T_p),s_tilde_nt_former(:,1:T_p),vs_tilde_nt_former(:,1:T_p),T_p,N_p,damping,tt,HbX_loc(:,1:T_p),Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_nt_former=s_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
 
        
           [Z_nt(:,1:T_p),V_nt(:,1:T_p),V_bar_nt(:,1:T_p),b_hat,x_hat(:,1:T_p),v_b,v_x(:,1:T_p),V1_se]=Forw_iter_SE(1,prior_X_mean(:,1:T_p),prior_X_var(:,1:T_p),prior_H,Sigma_x(:,1:T_p),R_x(:,1:T_p),Sigma_b,R_b,s_tilde_nt(:,1:T_p),X,T_p,Hr(1:K_p,:),K_p,V_nt_former(:,1:T_p),b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc(:,1:T_p),Hb_loc,HbX_loc(:,1:T_p),Hr_loc(1:K_p,:),V1_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T);
           V_nt_former=V_nt;
           b_hat_former=b_hat;
 
           [Z_kt(1:K_p,1:T_p),V_kt(1:K_p,1:T_p),V_bar_kt(1:K_p,1:T_p),q_hat(1:K_p,:),u_hat(:,1:T_p),v_q(1:K_p,:),v_u(:,1:T_p),V2_se]=Forw_iter_SE(2,Z_nt(:,1:T_p),V_nt(:,1:T_p),prior_Q,Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,1:T_p),X,T_p,Hr(1:K_p,:),K_p,V_kt_former(1:K_p,1:T_p),q_hat_former(1:K_p,:),N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc(:,1:T_p),Hb_loc,HbX_loc(:,1:T_p),Hr_loc(1:K_p,:),V2_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T); 
           V_kt_former=V_kt; 
           q_hat_former=q_hat;
 
        else
            [Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,:),vs_tilde_kt(1:K_p,:),v_tilde_2(1:K_p,:),vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(2,Rec_Y_pilot,Y_noise(1:K_p,:),Z_kt(1:K_p,:),V_kt(1:K_p,:),q_hat(1:K_p,:),u_hat,v_q(1:K_p,:),v_u,s_tilde_kt_former(1:K_p,:),vs_tilde_kt_former(1:K_p,:),T_p,N_p,damping,tt,HbX_loc,Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
 
            [Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,vs_tilde_nt,v_tilde_1,vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(1,R_u,Sigma_u,Z_nt,V_nt,b_hat,x_hat,v_b,v_x,s_tilde_nt_former,vs_tilde_nt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_nt_former=s_tilde_nt;
            vs_tilde_nt_former=vs_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
 
        
           [Z_nt,V_nt,V_bar_nt,b_hat,x_hat,v_b,v_x,V1_se]=Forw_iter_SE(1,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,X,T_p,Hr(1:K_p,:),K_p,V_nt_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc(1:K_p,:),V1_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T);
           V_nt_former=V_nt;
           b_hat_former=b_hat;
 
           [Z_kt(1:K_p,:),V_kt(1:K_p,:),V_bar_kt(1:K_p,:),q_hat(1:K_p,:),u_hat,v_q(1:K_p,:),v_u,V2_se]=Forw_iter_SE(2,Z_nt,V_nt,prior_Q,Sigma_u,R_u,Sigma_q(1:K_p,:),R_q(1:K_p,:),s_tilde_kt(1:K_p,:),X,T_p,Hr(1:K_p,:),K_p,V_kt_former(1:K_p,:),q_hat_former(1:K_p,:),N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc(1:K_p,:),V2_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T); 
           V_kt_former=V_kt; 
           q_hat_former=q_hat;
 
 
            diff(tt,1)=norm(x_dec-b_hat*x_hat,'fro')^2/norm(b_hat*x_hat,'fro')^2;
           if diff(tt,1)<1e-5
               break;
           end           
        end
        end %  outer loop stops (pilot parts) 
 
            
         
        A_eq=(b_hat*x_hat).'; % T*N
        prior_X_mean_new=zeros(N,K);
        prior_X_var_new=ones(N,K);
        prior_H_mean_new=zeros(T,N);
        prior_H_var_new=ones(T,N);
        [q_hat_new,v_q_new,u_hat_new,v_u_new]=BiGAMP(prior_X_mean_new,prior_X_var_new,prior_H_mean_new,prior_H_var_new,A_eq,T,damping_r,Rec_Y.',Y_noise.',Hr_loc.',var_channel,iter,Hr.',K_p);
        q_hat=q_hat_new.';
        v_q=v_q_new.';
        u_hat=u_hat_new';
        v_u=v_u_new';
%         


          
        for tt=1:4
             [Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt,vs_tilde_kt,v_tilde_2,vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(2,Rec_Y,Y_noise,Z_kt,V_kt,q_hat,u_hat,v_q,v_u,s_tilde_kt_former,vs_tilde_kt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_kt_former=s_tilde_kt;
            vs_tilde_kt_former=vs_tilde_kt;
   
            [Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,vs_tilde_nt,v_tilde_1,vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(1,R_u,Sigma_u,Z_nt,V_nt,b_hat,x_hat,v_b,v_x,s_tilde_nt_former,vs_tilde_nt_former,T_p,N_p,damping,tt,HbX_loc,Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T);
            s_tilde_nt_former=s_tilde_nt;
            vs_tilde_nt_former=vs_tilde_nt;
            x_hat_former=x_hat;
            x_dec=b_hat_former*x_hat_former; % the joint estimation
      
        
           [Z_nt,V_nt,V_bar_nt,b_hat,x_hat,v_b,v_x,V1_se]=Forw_iter_SE(1,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde_nt,X,T_p,q_hat_new.',K,V_nt_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,0.7,tt,S_loc,Hb_loc,HbX_loc,Hr_loc,V1_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T);
           V_nt_former=V_nt;
           b_hat_former=b_hat;
 
 
           [Z_kt,V_kt,V_bar_kt,q_hat,u_hat,v_q,v_u,V2_se]=Forw_iter_SE(2,Z_nt,V_nt,prior_Q,Sigma_u,R_u,Sigma_q,R_q,s_tilde_kt,X,T_p,q_hat_new.',K,V_kt_former,q_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,0.7,tt,S_loc,Hb_loc,HbX_loc,Hr_loc,V2_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T); 
           V_kt_former=V_kt; 
 
           q_hat_former=q_hat;
 
        end
        
 
   
 
         
        
       A_eq=(b_hat*x_hat).'; % T*N
        prior_X_mean_new=zeros(N,K);
        prior_X_var_new=ones(N,K);
        prior_H_mean_new=zeros(T,N);
        prior_H_var_new=ones(T,N);
        [q_hat_new,v_q_new,u_hat_new,v_u_new]=BiGAMP(prior_X_mean_new,prior_X_var_new,prior_H_mean_new,prior_H_var_new,A_eq,T,damping_r,Rec_Y.',Y_noise.',Hr_loc.',var_channel,iter,Hr.',K_p);
        q_hat=q_hat_new.';
        v_q=v_q_new.';
        u_hat=u_hat_new';
        v_u=v_u_new';
        

         chi_x=1;
         chi_b=1/M;
         chi_r=1/N;
         
         Sigma_u_se=sum(sum(Sigma_u))/sum(sum(Sigma_u~=0));
         chi_z2=N*M*chi_r*chi_b*chi_x;
         chi_z1=M*chi_b*chi_x;
         vars_x=0;
         vars_b=0;
         vars_r=0;
         V1=M*(chi_x*chi_b-vars_x*vars_b);
         vars_z1=chi_z1-V1+(V1)^2/(V1+Sigma_u_se);
         V2=N*(M*chi_b*chi_x*chi_r-vars_z1*chi_r);
         
         
         
          for tt=1:20
               % state evolution (two layers)
              vars_z2=chi_z2-V2+(V2)^2/(V2+var_noise);
              vars_z1=chi_z1-V1+(V1)^2/(V1+Sigma_u_se);
              
              Sigma_x_se=M^2*(chi_x*chi_b-vars_x*vars_b)^2/(N*vars_b*(vars_z1-M*vars_x*vars_b));
              Sigma_b_se=M^2*(chi_x*chi_b-vars_x*vars_b)^2/(T*vars_x*(vars_z1-M*vars_x*vars_b));
              Sigma_r_se=N^2*(M*chi_r*chi_x*chi_b-vars_z1*vars_r)^2/(T*vars_z1*(vars_z2-N*vars_z1*vars_r));
              Sigma_u_se=N^2*(M*chi_r*chi_x*chi_b-vars_z1*vars_r)^2/(T*vars_r*(vars_z2-N*vars_z1*vars_r));
              
              vars_x=(var_channel^2/(var_channel+Sigma_x_se));
              vars_b=(var_channel^2/(var_channel+Sigma_b_se));
              vars_r=(var_channel^2/(var_channel+Sigma_r_se));
              
              V1=M*(chi_x*chi_b-vars_x*vars_b);
              V2=N*(M*chi_b*chi_x*chi_r-vars_z1*chi_r);
          end
 
          
        
      se_x=(0.2*Sigma_x_se)/(0.2+Sigma_x_se);
      se_b=(0.2*Sigma_b_se)/(0.2+Sigma_b_se);
      se_r=(0.2*Sigma_r_se)/(0.2+Sigma_r_se);
 
         
         
end