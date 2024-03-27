function [out] = InterLeaver(cod )
%  InterLeaver is the location change
% modu_out is the modulation bits
% out is the location change sequence
N=length(cod);
loopMax=N/2;
s=fix(sqrt(N/2));
len=N;
genTime=0;

while (len~=0)
    len=N;
    genTime=genTime+1;
    loop=0;
    [~,a]=sort(rand(1,N));
    while (len~=0&&loop<loopMax)
%         [len, out]=mysort_new(N,s,a);
          stopflag=0;
          p=2;%p is the location to place the figure
          while(p<=N)
              cflag=1;
              next_flag=0;
              q=p;
              while(q<=N)
                  search=p-1;
                  while((p-search<=s)&&(search>=1))
                      cflag= ( abs(a(1,q) - a(1,search))>s )&& cflag;
                      search=search-1;
                  end
                  if(cflag~=1)
                      if(q==N)
                          stopflag=1;
                          break;
                      end
                      q=q+1;
                      next_flag=next_flag+1;
                      cflag=1;
                      continue
                  else
                      if( p == N && q==N )
                          len=0;
                          out=a;
                          return
                      end
                      if(next_flag~=0)
                          tmp=a(1,q);
                          a(1,q) = a(1,p) ;
                          a(1,p) = tmp ;
                      end
                      break;
                  end
              end
              if (stopflag ==1 )
                  len=N-p+1;
                  out=[a(1,p:N) a(1,1:p-1)] ;
                  break;
              end
              p=p+1;
          end
          
          a=out;
          loop=loop+1;
              
    end
    
    
end   




end

