function [post_mean,post_var]=VAMP(obser,measurement,pri_mean,pri_var,var_noise)

[length_row,length_col]=size(measurement);
% intialization
[U,S,V]=svd(measurement);
y_flag=1/var_noise*S'*U'*obser;

%% First method for initialization using LMMSE
gama_2(1,1)=1/pri_var;
x2_mean(:,1)=var_noise*inv(measurement'*measurement*var_noise+gama_2(1,1)*eye(length_col))*measurement'*obser;
x2_alpha=gama_2(1,1)/length_col*trace(inv(measurement'*measurement*var_noise+gama_2(1,1)*eye(length_col)));
x2_eta=gama_2(1,1)/x2_alpha;
gama_1(1,2)=x2_eta-gama_2(1,1);
r1(:,2)=x2_mean(:,1)*x2_eta/gama_1(1,2);
%%

% VAMP
% Start the iteration
T_NUM=15;
% damping=0.1; % The damping factor in the former iteration time
for t=2:1:T_NUM
    % Denoising
    % compute the belief of x1
    x1_eta(1,t)=(gama_1(1,t)+1/pri_var);
    x1_mean(:,t)= 1/x1_eta(1,t)*(r1(:,t)*gama_1(1,t)+pri_mean/pri_var);
%     if t~=2
%         x1_mean(:,t)=damping*x1_mean(:,t-1)+(1-damping)*x1_mean(:,t);
%     end
    % compute the message from x1 to delta
    gama_2(1,t)=x1_eta(1,t)-gama_1(1,t); 
    r2(:,t)=(x1_mean(:,t)*x1_eta(1,t)-r1(:,t)*gama_1(1,t))/gama_2(1,t);
    

    D=inv(1/var_noise*S'*S+gama_2(1,t)*eye(length_col));
    % LMMSE estimation
    % compute the belief of x2           
    x2_mean(:,t)=V*D*(y_flag+gama_2(1,t)*V'*r2(:,t));
    x2_eta_inv(1,t)=sum(diag(D)*gama_2(1,t))/length_row;
    r1(:,t+1)=(x2_mean(:,t)-x2_eta_inv(1,t)*r2(:,t))/(1-x2_eta_inv(1,t));
    gama_1(1,t+1)=gama_2(1,t)*(1-x2_eta_inv(1,t))/x2_eta_inv(1,t);
%     gama_1(1,t+1)=damping*gama_1(1,t)+(1-damping)*gama_1(1,t+1);
    
    i_vari(t,1)=0.0;
	a_fact(t,1)=0.0;
    i_vari(t,1)=norm(r1(:,t+1)-r1(:,t))^2;
    a_fact(t,1)=norm(r1(:,t))^2;
    if i_vari(t,1)<a_fact(t,1)*1e-12
        break;
    end
end
% End the iteration
%%
post_var=1/(1/pri_var+gama_1(1,t));
post_mean=post_var*(pri_mean/pri_var+r1(:,t)*gama_1(1,t));

end