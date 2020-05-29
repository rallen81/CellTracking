function y=mynormal(x,sigma,mu)

if sigma~=0
y=(1/(sqrt(2*pi)*sigma))*exp(-((x-mu).^2)/(2*sigma^2));
else
    if x==mu
        y=(1/(sqrt(2*pi)*sigma));
    end
end 


end

