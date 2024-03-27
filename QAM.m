function [ sym,qam_symbol ] = QAM(ori_bits )
% QAM modulation
P=2; % QAM bits
qam_symbol=zeros(2^P,1);
USED_SUB_NUM=length(ori_bits)/P;
sym=zeros(USED_SUB_NUM,1);
% Constellation points for 4-QAM
    for k=0:1
        for t=0:1
            qam_symbol(2*k+t+1)=(1-2*t)/sqrt(2)+1i*(1-2*k)/sqrt(2);
        end
    end
for i=0:USED_SUB_NUM-1
    code=ori_bits(P*i+1:P*(i+1),1);
    sym(i+1,1)=qam_symbol(2*code(1)+code(2)+1);
end
end

