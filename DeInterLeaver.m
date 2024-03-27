function deinterleavered_bits = DeInterLeaver( prob_bits,out )
% This function is the deinterleaver
% rec_awgny is the output ot the channel
% out is the location change of the interleaver

N=length(out);
for i=1:N
    deinterleavered_bits(out(1,i),:)=prob_bits(i,:);
end
end