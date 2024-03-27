function [Hs_mean,Hs_var]=ChannelEstimation(s_left,rhos_left,K,M,N,noise_var,damping,rec_y,h_prior_mean,h_prior_var)
% This function is to estimate channel using VAMP
% s_left and rho_left are the means and variances of the message
% transmitted from signal part to channel part

%% Compute the message from f_y^k to b_km
b_left=zeros(K,M);  
rhob_left=zeros(K,M);
sigma_g=sum(sum(rhob_left))/(M*K);

%%
% Initialization   
sma_num = 1e-6;
y_right=zeros(N,M);
gamma_w_right = sma_num;
gamma2 = sma_num*ones(M,1);
r2= zeros(K,M);
A=rand(K,N); 
%%
% Iteration starts
T_max=15;
stop=zeros(M,T_max);
for iter=1:T_max
    % Extrinsic part: message g_m to delta_s
    m=1:M;
    gamma_w_left=1/sigma_g-gamma_w_right;
    y_left(:,m)=(b_left(:,m)/sigma_g-gamma_w_right*y_right(:,m))/gamma_w_left;
    %%
    % LMMSE of Hs_2: the belief of Hs_2
    eta2=zeros(M,1);
    Hs_bar_2=zeros(N,M);
    for m=1:M 
        cov_hs2=inv(gamma_w_left*A'*A+gamma2(m,1)*eye(N));
        eta2(m,1)=N/trace(cov_hs2);
        Hs_bar_2(:,m)=cov_hs2*(gamma_w_left*A'*y_left(m,1)+gamma2(m,1)*r2(:,m));
    end
    %% the message from h^s_2 to delta_h
    gamma1=eta2-gamma2;
    r1=(repmat(eta2.',N,1).*Hs_bar_2-repmat(gamma2.',N,1).*r2)./(repmat(gamma1.',N,1));

    %% the belief of h1
    eta1=gamma1+1./h_prior_var;
    Hs_bar_1=(repmat(gamma1.',N,1).*r1+repmat((h_prior_mean/h_prior_var).',N,1));

    %% the message from h^s_1 to delta_h
    gamma2=eta1-gamma1;
    r2=(repmat(eta1.',N,1).*Hs_bar_1-repmat(gamma1.',N,1).*r1)./(repmat(gamma2.',N,1));

    %% the message to MMSE modeule: the belief of g_m
    for m=1:M 
        cov_z=A*inv(gamma2(m,1)*eye(N)+gamma_w_left*A'*A)*A';
        eta_z(m,1)=trace(cov_z)/N;
        G_bar(:,m)=A*inv(gamma2(m,1)*eye(N)+gamma_w_left*A'*A)*(gamma2(m,1)*r2(:,m)+gamma_w_left*A'*y_left(:,m));
    end
    eta=sum(eta_z)/M;

    %% the message from delta_s to g_m
    gamma_w_right=1/eta-gamma_w_left;
    for m=1:M
        y_right(:,m)=(G_bar(:,m)/eta-gamma_w_left*y_left(:,m))/gamma_w_right;
    end
    %% Iteration stop
    stop(:,iter)=r1;
    if iter~=1
        i_vari(iter,1)=0.0;
        a_fact(iter,1)=0.0;
        i_vari(iter,1)=norm(stop(:,iter)-stop(:,iter-1))^2;
        a_fact(iter,1)=norm(stop(:,iter-1))^2;
        if i_vari(iter,1)<a_fact(iter,1)*1e-12
            break;
        end
    end
end
%%
post_var=1/(1/h_prior_var+gamma1);
post_mean=post_var*(h_prior_mean/h_prior_var+r1*gamma1);

end