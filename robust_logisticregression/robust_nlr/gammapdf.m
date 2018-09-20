function y = gammapdf(x,k,theta)

x = abs(x);

y = ((x.^(k-1)).*exp(-x/theta))/(theta^k * factorial(k-1));
