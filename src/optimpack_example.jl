using OptimPackNextGen

function chi2_fg(x, g, data) # criterion function plus its gradient w/r x
  f = sum((x-data).^2)
  g[:] = 2*(x-data)
  return f
end

data = rand(100).^2;
x_start = zeros(100);

crit = (x,g)->chi2_fg(x, g, data);
x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80);
