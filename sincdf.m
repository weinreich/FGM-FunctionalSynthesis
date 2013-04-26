function cdf = sincdf(x,n)
% Implement cdf for pdf(x,n) = Zsin[x]^(n-2).
% Got result from Mathematica!

if  (((n == 2) || (n == 3)) && (x == 90))
    %some kind of matlab issue causes hypergeom() to blow up when last
    %argument is zero and n is 2 or 3. Happily the function equals 1 at
    %that point.
     cdf = gamma(n/2)*(sqrt(pi)*gamma((n-1)/2)/(2*gamma(n/2))-cosd(x))/(sqrt(pi)*gamma((n-1)/2));
else
    cdf = gamma(n/2)*(sqrt(pi)*gamma((n-1)/2)/(2*gamma(n/2))-cosd(x)*hypergeom([1/2 3/2-n/2],3/2,cosd(x)^2))/(sqrt(pi)*gamma((n-1)/2));
end
