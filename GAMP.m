function [post_mean,post_var]=GAMP(obser,measurement,pri_mean,pri_var,var_noise)


[length_row,length_col]=size(measurement);


% GAMP
% intialization
Z_a(:,1)=obser;
V_a(:,1)=ones(length_row,1);
for k=1:length_col
    h_mean(k,2)=pri_mean;
    h_var(k,2)=pri_var;
end

% Start the iteration
T_NUM=15;
damping=1; % The damping factor in the former iteration time
for t=2:1:T_NUM
    for j=1:length_row
        temp1=0;
        temp2=0;
        for k=1:length_col
            temp2=temp2+norm(measurement(j,k))^2*h_var(k,t);
            temp1=temp1+measurement(j,k)*h_mean(k,t);
        end
%         if t~=2
%             V_a(j,t)=damping*V_a(j,t-1)+(1-damping)*temp2;
%         else
            V_a(j,t)=temp2;
%         end
        Z_a(j,t)=temp1-V_a(j,t)/(var_noise+V_a(j,t-1))*(obser(j)-Z_a(j,t-1));
    end
    
    for k=1:length_col
        temp2=0;
        temp3=0;
        for j=1:length_row
            temp2=temp2+norm(measurement(j,k))^2/(var_noise+V_a(j,t));
            temp3=temp3+conj(measurement(j,k))*(obser(j)-Z_a(j,t))/(var_noise+V_a(j,t));
        end
        S_j(k,t)=1/temp2;
%         if t~=2
%             h_mean_tilde(k,t)=damping*h_mean_tilde(k,t-1)+(1-damping)*h_mean(k,t);
%         else
%             h_mean_tilde(k,t)=h_mean(k,t);
%         end
        R_j(k,t)=h_mean(k,t)+S_j(k,t)*temp3;
    end
    
    % Denoising step
    for k=1:length_col
        h_var(k,t+1)=1/(1/h_var(k,t)+1/S_j(k,t));
        h_mean(k,t+1)=R_j(k,t)/S_j(k,t)*h_var(k,t+1);
%         h_mean(k,t+1)=damping*h_mean(k,t)+(1-damping)*h_mean(k,t+1);
    end
    
    i_vari(t,1)=0.0;
	a_fact(t,1)=0.0;
    i_vari(t,1)=norm(h_mean(:,t+1)-h_mean(:,t))^2;
    a_fact(t,1)=norm(h_mean(:,t))^2;
    if i_vari(t,1)<a_fact(t,1)*1e-12
        break;
    end
end

post_var=h_var(:,t);
post_mean=h_mean(:,t);

end