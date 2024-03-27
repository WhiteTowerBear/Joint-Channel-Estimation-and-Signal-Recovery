function [Sigma_u,R_u,Sigma_q,R_q,s_tilde,vs_tilde,v_tilde,vars_x,vars_b,vars_r,vars_z1]=Back_iter_SE(flag_back,Rec_Y,var_noise,Z,V,q_hat,u_hat,v_q,v_u,s_tilde_former,vs_tilde_former,T_p,N_p,damping,tt,HbX_loc,Hb,X,V1_se,V2_se,vars_z1,vars_x,vars_b,vars_r,vs1_se_former,vs2_se_former,M,N,K,T)
% Code for paper "Message passing for joint CE and data detection in RIS-based MISO system"
% Code written by Wei Li
% Jan 14, 2021

% X0, H0 and Q0 are intialization of estimated X, H and Q
% prior_X, prior_H and prior_Q are priors of estimated X, H and Q (The 
%    first element is prior mean, and the second element is prior variance)
% Rec_Y is the received signal 
% var_noise is the variance of AWGN 
% Z is the product of H^b and X (or U and Q), V is variance of the product

v_tilde=1 ./(1./V+1./var_noise);
 
    if flag_back==2
        Z2_se=sum(sum(Z))/sum(sum(Z~=0));
        vars_z2=Z2_se^2;
        chi_z2=V2_se+vars_z2;
        v2_tilde_se=chi_z2-vars_z2; % state evolution
        vs2_se=(V2_se-v2_tilde_se)/(V2_se^2);
        if tt~=1
            vs2_se=damping*vs2_se+(1-damping)*vs2_se_former;
        end
         Sigma_r_se=1/(T*vars_z1*vs2_se);
         Sigma_u_se=1/(K*vars_r*vs2_se);
         vars_r=(1/(1+Sigma_r_se));
         vars_z1=(1/(1+Sigma_u_se));
    else
        Z1_se=sum(sum(Z))/sum(sum(Z~=0));
        vars_z1=Z1_se^2;
        chi_z1=V1_se+vars_z1;
        v1_tilde_se=chi_z1-vars_z1;
        vs1_se=(V1_se-v1_tilde_se)/(V1_se^2);
        if tt~=1
            vs1_se=damping*vs1_se+(1-damping)*vs1_se_former;
        end
         Sigma_b_se=1/(T*vars_x*vs1_se);
         Sigma_x_se=1/(N*vars_b*vs1_se);
         vars_x=(1/(1+Sigma_x_se));
         vars_b=(1/(1+Sigma_b_se));
    end
 

z_tilde=v_tilde.*(Z./V+Rec_Y./var_noise);

 
s_tilde= ( z_tilde - Z)./V;
if tt~=1
    s_tilde=damping*s_tilde+(1-damping)*s_tilde_former;
end
vs_tilde=(V - v_tilde) ./ (V.^2);
 
if tt~=1
    vs_tilde=damping*vs_tilde+(1-damping)*vs_tilde_former;
end

 v_tilde(find(isnan(v_tilde)==1))=0; 
z_tilde(find(isnan(z_tilde)==1))=0;
s_tilde(find(isnan(s_tilde)==1))=0; 
vs_tilde(find(isnan(vs_tilde)==1))=0;

Sigma_u=1./(abs(q_hat.').^2*vs_tilde);
% Sigma_u_scalar=1./(sum(sum(abs(q_hat).^2))*vs_tilde_scalar/size(q_hat,1));
R_u=u_hat.*(1-Sigma_u.*(v_q.'*vs_tilde))+Sigma_u.*(q_hat'*s_tilde);
if flag_back==2
    Sigma_u=Sigma_u.*HbX_loc;
    R_u=R_u.*HbX_loc;
    u_hat(1:N_p,1:T_p)=Hb(1:N_p,:)*X(:,1:T_p);
    v_u(1:N_p,1:T_p)=0;
    u_hat=u_hat.*HbX_loc;
    v_u=v_u.*HbX_loc;
end
Sigma_u(find(isnan(Sigma_u)==1))=0; 
R_u(find(isnan(R_u)==1))=0;

Sigma_q=1./(vs_tilde* abs(u_hat.').^2);
R_q=q_hat.*(1-Sigma_q.*(vs_tilde*v_u.'))+Sigma_q.*(s_tilde*u_hat');



end