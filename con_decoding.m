function [detect_c]= con_decoding( llr_ext )
% con_decoding aims at decoding convolutional codes
% prior_prob is the prior probablity of each received symbols
% calculate app_prob: posterior probability
% calculate 
% [llr_app,detect_u] 
% N_length=size(rec_prob);
% llr_capp stores the llr probability of coded bits

% trellis grap

% trellis=zeros(64,5); trellis(:,1)is the former state, trellis(:,2)is the
% next state corresponding to the input 0, trellis(:,3) is the next state
% corresponding to the input 1, trellis(:,4) is the output corresponding to
% the input 0, trellis(:,5) is the output corresponding to the input 1.
prob_all(:,2)=1./(1+exp(-llr_ext));
prob_all(:,1)=1-prob_all(:,1);
llr_c1=prob_all;
input_num=length(llr_c1);
trellis = poly2trellis(7,[133 171 165]);
LLR=BCJR_conv(llr_c1,trellis);
  for m=1:length(LLR)-6    
      if LLR(m)<0
          detect_c(m,1)=0;
      else
          detect_c(m,1)=1;
      end
   end 

end


