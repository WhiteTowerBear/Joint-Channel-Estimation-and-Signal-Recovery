function [Z,V,V_bar,b_hat,x_hat,v_b,v_x]=Single_Forw_iter(prior_X_mean,prior_X_var,prior_H_mean,prior_H_var,Sigma_x,R_x,Sigma_b,R_b,s_tilde,Hr,K_p,V_former,b_hat_former,damping,tt,S_loc,X,T_p)
% Code for paper "Message passing for joint CE and data detection in RIS-based MISO system"
% Code written by Wei Li
% Jan 14, 2021

% prior_X, prior_H and prior_Q are priors of estimated X, H and Q (The 
%    first element is prior mean, and the second element is prior variance)
% Rec_Y is the received signal 
% var_noise is the variance of AWGN 
% Z is the product of H^b and X (or U and Q), V is variance of the product


v_x=1./(1./Sigma_x+1./prior_X_var);
v_x=v_x.*S_loc;
x_hat=v_x.*(R_x./Sigma_x+prior_X_mean./prior_X_var);
x_hat=x_hat.*S_loc;
x_hat(:,1:T_p)=X(:,1:T_p);
v_x(:,1:T_p)=0;
v_b=1./(1./Sigma_b+1./prior_H_var);
b_hat=v_b.*(R_b./Sigma_b+prior_H_mean./prior_H_var);



    if tt~=1
       b_hat=damping*b_hat+(1-damping)*b_hat_former;
    end
    if K_p~=0
        b_hat(1:K_p,:)=Hr(1:K_p,:); % The known channel part
        v_b(1:K_p,:)=0;
    end


V_bar=v_b*abs(x_hat).^2+abs(b_hat).^2*v_x;
% if tt~=1 
%     V_bar=damping*V_bar+(1-damping)*V_bar_former;
% end
Z_bar=b_hat*x_hat;
V=V_bar+v_b*v_x;
if tt~=1
    V=damping*V+(1-damping)*V_former;
end
Z=Z_bar-s_tilde.*V_bar; 


% if tt~=1
%     Z=damping*Z+(1-damping)*Z_former;
% end

end