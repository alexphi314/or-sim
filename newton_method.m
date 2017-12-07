function [root, counter] = newton_method(fh,f_prime, init,tolerance)
%Input an anonymous function handle, the derivative of that function, an initial guess, and the tolerance to implement the
%Newton method

%Defining initial conditions and performing the first iteration
fx = fh(init);
prime = f_prime(init);
ratio = fx./prime;
x = init;
counter = 1;

while abs(ratio) > tolerance
    x = x - ratio;
    fx = fh(x);
    prime = f_prime(x);
    ratio = fx./prime;
    %display(ratio);
    counter = counter + 1;
end

root = x;
end