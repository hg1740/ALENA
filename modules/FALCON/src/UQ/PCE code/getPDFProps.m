function [mean,std,skewness,kurtosis] = getPDFProps(x,fx)

mean     = trapz(x, x         .*fx);
variance = trapz(x,(x-mean).^2.*fx);
std      = sqrt(variance);
skewness = trapz(x,(x-mean).^3.*fx)/std^3;
kurtosis = trapz(x,(x-mean).^4.*fx)/std^4;