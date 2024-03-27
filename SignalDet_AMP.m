function [ llr_amp,estimated_codes,x_esti] = SignalDet_AMP( rec_y,H_eq,qam_symbol,M,var_noise)
% Input: the equivalent channel matrix, e.g., H_eq, the measurement
% Out put: the llr value of each coded bit

% constellation
num_conste=4; % QAM modulation 

% intialization
pri_mean=0;
pri_var=sum(norm(qam_symbol).^2)/num_conste;
% VAMP for data detection
[post_mean,post_var]=VAMP(rec_y,H_eq,pri_mean,pri_var,var_noise);

for i=1:M
    for j=1:num_conste
        p_com(i,j)=exp(-(norm(qam_symbol(j,1)-post_mean(i,1))^2)/post_var);
    end
end
p_sum=(sum(p_com,2));
p_com=p_com./p_sum;

%%
[~,index]=max(p_com,[],2);
kk=1:1:length(p_com);
x_esti(kk,1)=qam_symbol(index(kk));
%%
p_com_re=kron(p_com,[1;1]);
prob_bits=zeros(M*num_conste/2,2); % stores the probability of each coded bit
prob_bits(1:2:end,1)=p_com_re(1:2:end,1)+p_com_re(1:2:end,2); % the bit in odd location 
prob_bits(1:2:end,2)=p_com_re(1:2:end,3)+p_com_re(1:2:end,4); 
prob_bits(2:2:end,1)=p_com_re(2:2:end,1)+p_com_re(2:2:end,3); % the bit in even location 
prob_bits(2:2:end,2)=p_com_re(2:2:end,2)+p_com_re(2:2:end,4); 

llr_amp=log(prob_bits(:,1)./prob_bits(:,2));
 
estimated_codes=(llr_amp<0);
end
