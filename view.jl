using PyPlot

function imview(x;title="")
    fig = figure(title,figsize=(10,10))
    imshow(x, cmap=ColorMap("coolwarm"), interpolation="none");
    tight_layout()
end

function imsurf(z, title="")
    fig = figure(title,figsize=(10,10))
    x = collect(1:size(z,1))
    y = collect(1:size(z,2))
    zz = copy(z)
    zz[isnan.(z)]=0
    surf(x,y,z, facecolors=get_cmap("coolwarm")(zz/maximum(zz)), linewidth=0.25, rstride=4, edgecolors="k", cstride=4,antialiased=false, shade=false)
    tight_layout()
end
