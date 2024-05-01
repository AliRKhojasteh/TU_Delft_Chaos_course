function xs = Runge_slave(xs,xm)
   global h ga Gam Om c
   
   x  = zeros(size(xs));
   c1 = zeros(size(xs));
   c2 = zeros(size(xs));
   c3 = zeros(size(xs));
   c4 = zeros(size(xs));

   n = length(xs);
   for i = 1:n; x(i) = xs(i); end
   f = equations2(x,xm);
   for i = 1:n; c1(i) = h*f(i); end

   for i = 1:n; x(i) = xs(i) + c1(i)/2; end
   f = equations2(x,xm);
   for i = 1:n; c2(i) = h*f(i); end

   for i = 1:n;  x(i) = xs(i) + c2(i)/2; end
   f = equations2(x,xm);
   for i = 1:n;  c3(i) = h*f(i); end

   for i = 1:n;  x(i) = xs(i) + c3(i); end
   f = equations2(x,xm);
   for i = 1:n;  c4(i) = h*f(i); end
   
   for i = 1:n
       xs(i) = xs(i) + (c1(i) + 2*c2(i) + 2*c3(i) + c4(i))/6;
   end
end   


function f = equations2(xs,xm)
    global ga Gam Om c
    f = zeros(3,1);
    f(1) = xs(2);
    f(2) = -ga * xs(2) - sin(xs(1)) + Gam*cos(Om*xs(3)) + c * (sin(xs(1)) - sin(xm(1))); 
    f(3) = 1;
end
