using PyPlot

function imview(x;title="",zoom=1, color="gray")
    if zoom < 1
        zoom = 1
    end
    fig = figure(title,figsize=(10,10))
    clf()
    nx = size(x,1)
    ny = size(x,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))
    imshow(x[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y], cmap=ColorMap(color), interpolation="none");
    tight_layout()
end

function imsurf(z; title="", zoom=1, color="coolwarm")
    if zoom < 1
        zoom = 1
    end
    fig = figure(title,figsize=(10,10))
    clf()
    axis("equal")
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
    surf(x,y,zoomed, facecolors=get_cmap(color)(zoomed/maximum(zoomed)), linewidth=0.25, rstride=4, edgecolors="k", cstride=4,antialiased=false, shade=false)
    tight_layout()
end

function imview_add(x;zoom=1, color="gray") #will overplot whatever you're drawing
    if zoom < 1
        zoom = 1
    end
    nx = size(x,1)
    ny = size(x,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))
    imshow(x[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y], cmap=ColorMap(color), interpolation="none");
    tight_layout()
end
