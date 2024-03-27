function [ tran_x,pilot_location, data_location] = pilot( M,D_p,M_p, sym,pilot_symbol )
% M_p is the number of pilots
tran_x=zeros(M+M_p,1);

gg=0:1:(M_p-1);
tran_x((D_p+1)*gg+1,1)=pilot_symbol;
for kk=0:1:(M_p-2)
    tran_x((D_p+1)*kk+2:(D_p+1)*kk+D_p+1,1)=sym(D_p*kk+1:D_p*kk+D_p,1);
end
pilot_location=zeros(M_p,1);
pilot_location(gg+1,1)=(D_p+1)*gg+1;
data_location=zeros(M,1);
for kk=0:1:(M_p-2)
    data_location(D_p*kk+1:D_p*kk+D_p,1)=(D_p+1)*kk+2:(D_p+1)*kk+D_p+1;
end


end

