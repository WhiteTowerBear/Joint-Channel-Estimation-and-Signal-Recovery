function [x_hat,a_hat]=BiGAMP_lite(Y,M,L,N,var_noise,x_prior,vx_prior,a_prior,va_prior,damping)
% [M,L]=size(Y);
% [M,N]=size(A);
% [N,L]=size(X);
v_hat(:,:,1)=zeros(M,L);
v_x=1;
x_hat=(sqrt(vx_prior)*(randn(N,L)));
v_a=1;
a_hat=sqrt(va_prior)*(randn(M,N));
T_max=15;
for t=2:T_max
    Ga=N/(norm(a_hat,'fro')^2);
    Gx=N/(norm(x_hat,'fro')^2);
    u_hat(:,:,t)=Y-a_hat*x_hat;
    v_p_bar=(v_x/(M*Ga)+v_a/(L*Gx))*N;
    v_p=v_p_bar+N*v_a*v_x;
    v_hat(:,:,t)=u_hat(:,:,t)+v_p_bar/(v_p+var_noise)*v_hat(:,:,t-1);
    v_r=Ga*(v_p+var_noise);
    v_q=Gx*(v_p+var_noise);
    v_x=1/(1/v_r+1/vx_prior);
    x_hat=v_x/v_r*((1-M*v_a*Ga)*x_hat+Ga*a_hat'*v_hat(:,:,t));
    v_a=1/(1/v_q+1/va_prior);
    a_hat=v_a/v_q*((1-L*v_x*Gx)*a_hat+Gx*v_hat(:,:,t)*x_hat');
   
    
    if t~=2
    i_vari(t,1)=norm(u_hat(:,:,t)-u_hat(:,:,t-1),'fro')^2;
    a_vari(t,1)=1e-5*norm(u_hat(:,:,t),'fro')^2;
        if i_vari<a_vari
            break;
        end
    end
end


end
