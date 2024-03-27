function [location,cor_matrix]=correlated_matrix(input_matrix)
% This function is t ogenerate a correlated sparse matrix
[N,M]=size(input_matrix);
Hb_loc=binornd(1,0.2,N,M); 
Hb=input_matrix.*Hb_loc;

cor_matrix=zeros(N,M);
Hb_part=triu(Hb,2)+tril(Hb,-2);
Hb_diag=Hb-Hb_part;
cor_matrix(1:min(N,M),1:min(N,M))=tril(Hb_diag(1:min(N,M),1:min(N,M)),-1)+tril(Hb_diag(1:min(N,M),1:min(N,M)),-1).'+diag(diag(Hb));
cor_matrix=cor_matrix+Hb_part;

location=zeros(N,M);
Hb_loc_part=triu(Hb_loc,2)+tril(Hb_loc,-2);
Hb_loc_diag=Hb_loc-Hb_loc_part;
location(1:min(N,M),1:min(N,M))=tril(Hb_loc_diag(1:min(N,M),1:min(N,M)),-1)+tril(Hb_loc_diag(1:min(N,M),1:min(N,M)),-1).'+diag(diag(Hb_loc));
location=location+Hb_loc_part;
end