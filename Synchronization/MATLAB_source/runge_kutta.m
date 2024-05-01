function xm = Runge(xm)
   global h ga Gam Om
   
   x  = zeros(size(xm));
   c1 = zeros(size(xm));
   c2 = zeros(size(xm));
   c3 = zeros(size(xm));
   c4 = zeros(size(xm));

   n = length(xm);
   for i = 1:n; x(i) = xm(i); end
   f = equations1(x);
   for i = 1:n; c1(i) = h*f(i); end

   for i = 1:n; x(i) = xm(i) + c1(i)/2; end
   f = equations1(x);
   for i = 1:n; c2(i) = h*f(i); end

   for i = 1:n;  x(i) = xm(i) + c2(i)/2; end
   f = equations1(x);
   for i = 1:n;  c3(i) = h*f(i); end

   for i = 1:n;  x(i) = xm(i) + c3(i); end
   f = equations1(x);
   for i = 1:n;  c4(i) = h*f(i); end
   
   for i = 1:n
       xm(i) = xm(i) + (c1(i) + 2*c2(i) + 2*c3(i) + c4(i))/6;
   end
end   

function f = equations1(x)
    global ga Gam Om
    f = zeros(3,1);
    f(1) = x(2);
    f(2) = -ga * x(2) - sin(x(1)) + Gam*cos(Om*x(3)); 
    f(3) = 1;
end

