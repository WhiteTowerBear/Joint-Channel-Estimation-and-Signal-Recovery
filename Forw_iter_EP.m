function [Z,V,V_bar,b_hat,x_hat,v_b,v_x]=Forw_iter_EP(flag_forw,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde,X,T_p,Hr,K_p,V_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc)
% Code for paper "Message passing for joint CE and data detection in RIS-based MISO system"
% Code written by Wei Li
% Jan 14, 2021

% prior_X, prior_H and prior_Q are priors of estimated X, H and Q (The 
%    first element is prior mean, and the second element is prior variance)
% Rec_Y is the received signal 
% var_noise is the variance of AWGN 
% Z is the product of H^b and X (or U and Q), V is variance of the product
 
if flag_forw==1
   [M,T]=size(Sigma_x);
   px=zeros(M,T,2); % 2 is the constellation points, BPSK
   px_ln(:,:,1)=(-(-1/sqrt(2)-R_x).^2./Sigma_x);
   px_ln(:,:,2)=(-(1/sqrt(2)-R_x).^2./Sigma_x);
   px1=1./(1+exp(px_ln(:,:,2)-px_ln(:,:,1)));
   px2=1-px1;
   joint_X_mean=1/sqrt(2)*px2-1/sqrt(2)*px1;
   joint_X_var=1/2*px2+1/2*px1-joint_X_mean.^2;
   joint_X_var=abs(joint_X_var);
   prior_X_var=1./(1./joint_X_var-1./Sigma_x);
   prior_X_var=abs(prior_X_var);
   prior_X_mean=prior_X_var.*(joint_X_mean./joint_X_var-R_x./Sigma_x);
 
end

v_x=1./(1./Sigma_x+1./prior_X_var);
x_hat=v_x.*(R_x./Sigma_x+prior_X_mean./prior_X_var);


if T_p~=0
    if flag_forw==2
        x_hat(find(isnan(x_hat)==1))=0; 
        v_x(find(isnan(v_x)==1))=0;
    end
end
 


v_b=1./(1./Sigma_b+1./prior_H(2));
b_hat=v_b.*(R_b./Sigma_b+prior_H(1)./prior_H(2));

if T_p~=0
    if flag_forw==1
        x_hat(:,1:T_p)=X(:,1:T_p); % The pilot data part
        v_x(:,1:T_p)=0;
        x_hat=x_hat.*S_loc;
        v_x=v_x.*S_loc;
        b_hat(1:N_p,:)=Hb(1:N_p,:);
        v_b(1:N_p,:)=0;
        b_hat=b_hat.*Hb_loc;
        v_b=v_b.*Hb_loc;
    else
        x_hat(1:N_p,1:T_p)=Hb(1:N_p,:)*X(:,1:T_p);
        v_x(1:N_p,1:T_p)=0;
        x_hat=x_hat.*HbX_loc;
        v_x=v_x.*HbX_loc;
    end
end


 

if flag_forw==2
    if tt~=1
       b_hat=damping*b_hat+(1-damping)*b_hat_former;
    end
    if K_p~=0
        b_hat(1:K_p,:)=Hr(1:K_p,:); % The known channel part
        b_hat=b_hat.*Hr_loc;
        v_b=v_b.*Hr_loc;
        v_b(1:K_p,:)=0;
    end
end

V_bar=v_b*abs(x_hat).^2+abs(b_hat).^2*v_x;
 
Z_bar=b_hat*x_hat;
V=V_bar+v_b*v_x;
if tt~=1
    V=damping*V+(1-damping)*V_former;
end
Z=Z_bar-s_tilde.*V_bar; 


Z(find(isnan(Z)==1))=0; 
V(find(isnan(V)==1))=0;

end