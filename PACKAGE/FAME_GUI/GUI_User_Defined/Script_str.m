% ¦r¦êÂàoutput
n_temp=length(temp);
for k=1:n_temp
    ctp=char(temp);
    eval( [ctp(k,:),';']);
end