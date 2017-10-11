function cog(x) #note: this takes a 2D array
xvals=[i for i=1:size(x,1)]
return vec([sum(xvals'*x) sum(x*xvals)]/sum(x))
end

function posmax(x) #note: this takes a 2D array
pos = ind2sub(size(x),indmax(x))
return [pos[1],pos[2]]
end

function recenter(x, newcenter) #note: this takes a 2D array
Δ=[(size(x,1)+1)/2, (size(x,2)+1)/2]-newcenter
return circshift(x,round.(Δ))
end
