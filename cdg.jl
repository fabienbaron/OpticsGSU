function cdg(x) #note: this takes a 2D array
xvals=[i for i=1:size(x,1)]
return [sum(xvals'*x) sum(x*xvals)]/sum(x)
end
