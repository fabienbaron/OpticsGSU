function gaussian2d(n,m,sigma)
g2d = [exp(-((X-(div(m,2)+1)).^2+(Y-(div(n,2)+1)).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
return g2d
end
