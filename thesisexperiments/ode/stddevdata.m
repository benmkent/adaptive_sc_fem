stddevfn = @(t,E) exp(-E *t) .* sqrt(1 - (sinc(t/pi).^2))

t = logspace(-2,2,1e3);
results = [0;t(:)];

for exponent = -1:1;
    s = stddevfn(t,10^(exponent));
    results = [results,[ exponent; s(:)]];
end

writematrix(results,'stddev-ode.dat');