function f = vaccination(Vaccination,TA,t,NP)
f = zeros(NP,length(t));
for jj = 1:NP
f(jj,:) = interp1(TA,Vaccination(jj,:),t);
end