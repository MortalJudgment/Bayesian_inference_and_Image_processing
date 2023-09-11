N = 10^6;
kmin= 1;    
kmax = 30;
posteriors = [];
for k = kmin:kmax
    posteriors = [posteriors 2^k/(2^k + N - 1)];
end
plot( (kmin:kmax), posteriors, '-','LineWidth',1);
hold on
plot( (kmin:kmax), posteriors, 'o');
xlabel('number of flips')
ylabel('posterior probability')

