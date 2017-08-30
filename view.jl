using PyPlot

function imview(x;title="",zoom=1)
    fig = figure(title,figsize=(10,10))
    clf()
    nx = size(x,1)
    ny = size(x,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))
    imshow(x[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y], cmap=ColorMap("coolwarm"), interpolation="none");
    tight_layout()
end

function imsurf(z; title="", zoom=1)
    fig = figure(title,figsize=(10,10))
    clf()
    nx = size(z,1)
    ny = size(z,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))
    zoomed = z[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y]
    zoomed[isnan.(zoomed)]=0
    x = collect(1:size(zoomed,1))
    y = collect(1:size(zoomed,2))
    surf(x,y,zoomed, facecolors=get_cmap("coolwarm")(zoomed/maximum(zoomed)), linewidth=0.25, rstride=1, edgecolors="k", cstride=1,antialiased=false, shade=false)
    tight_layout()
end
