function cod=codec(m)
% Convolutional coding
g1=[1 0 1 1 0 1 1];
g2=[1 1 1 1 0 0 1];
g3=[1 1 1 0 1 0 1];
m1=conv(m,g1);   
m2=conv(m,g2); 
m3=conv(m,g3);
l=length(m1); 
cod=zeros(3*l,1);
for i=1:l  
    cod(3*i-2)=m1(i);
    cod(3*i-1)=m2(i);
    cod(3*i)=m3(i);
end

for i=1:(3*l)
    if(mod(cod(i),2))==0
        cod(i)=0;
    else
        cod(i)=1;
    end
end

end  

