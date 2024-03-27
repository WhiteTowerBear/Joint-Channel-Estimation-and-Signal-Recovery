function [Z,V,V_bar,b_hat,x_hat,v_b,v_x,V_se]=Forw_iter_SE(flag_forw,prior_X_mean,prior_X_var,prior_H,Sigma_x,R_x,Sigma_b,R_b,s_tilde,X,T_p,Hr,K_p,V_former,b_hat_former,N_p,Hb,Z_nt_pilot,V_nt_pilot,damping,tt,S_loc,Hb_loc,HbX_loc,Hr_loc,V_se_former,chi_x,chi_b,chi_r,vars_x,vars_b,vars_r,vars_z1,M,N,K,T)
% Code for paper "Message passing for joint CE and data detection in RIS-based MISO system"
% Code written by Wei Li
% Jan 14, 2021

% prior_X, prior_H and prior_Q are priors of estimated X, H and Q (The 
%    first element is prior mean, and the second element is prior variance)
% Rec_Y is the received signal 
% var_noise is the variance of AWGN 
% Z is the product of H^b and X (or U and Q), V is variance of the product
 
v_x=1./(1./Sigma_x+1./prior_X_var);
% v_x_scalar=sum(sum(v_x))/sum(sum(v_x~=0));
x_hat=v_x.*(R_x./Sigma_x+prior_X_mean./prior_X_var);
if T_p~=0
    if flag_forw==2
        x_hat(find(isnan(x_hat)==1))=0; 
        v_x(find(isnan(v_x)==1))=0;
    end
end
 


v_b=1./(1./Sigma_b+1./prior_H(2));
% v_b_scalar=sum(sum(v_b))/sum(sum(v_b~=0));
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
 
if flag_forw==1
    V_se=M*(chi_x*chi_b-vars_x*vars_b);
else
    V_se=N*(M*chi_b*chi_x*chi_r-vars_z1*vars_r);
end
 
if tt~=1
    V=damping*V+(1-damping)*V_former;
    V_se=damping*V_se+(1-damping)*V_se_former;
end
Z=Z_bar-s_tilde.*V_bar; 


Z(find(isnan(Z)==1))=0; 
V(find(isnan(V)==1))=0;

end