function [Sigma_u,R_u,Sigma_q,R_q,s_tilde,vs_tilde]=Single_Back_iter(Rec_Y,var_noise,Z,V,q_hat,u_hat,v_q,v_u,s_tilde_former,vs_tilde_former,damping,tt)
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

z_tilde=v_tilde.*(Z./V+Rec_Y./var_noise);
 
s_tilde= ( z_tilde - Z)./V;
if tt~=1
    s_tilde=damping*s_tilde+(1-damping)*s_tilde_former;
end
vs_tilde=(V - v_tilde) ./ (V.^2);
if tt~=1
    vs_tilde=damping*vs_tilde+(1-damping)*vs_tilde_former;
end
 

Sigma_u=1./(abs(q_hat.').^2*vs_tilde);
Sigma_q=1./(vs_tilde* abs(u_hat.').^2);
R_u=u_hat.*(1-Sigma_u.*(v_q.'*vs_tilde))+Sigma_u.*(q_hat'*s_tilde);
R_q=q_hat.*(1-Sigma_q.*(vs_tilde*v_u.'))+Sigma_q.*(s_tilde*u_hat');

end