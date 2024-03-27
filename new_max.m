function c=new_max(a,b)
% c=max(a,b)+log(1+exp(-abs(a-b)));


% if a==1.000000000000000e-100
%     c=b;
% elseif b==1.000000000000000e-100
%     c=a;
% else
%         c=log(exp(a)+exp(b));
%  end


c=max(a,b);
c=c+log(1+exp(-abs(a-b)));
end