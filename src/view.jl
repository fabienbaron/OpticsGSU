using PyPlot

function imview(x;title="",zoom=1, color="gray", newfig = true)
    if zoom < 1
        zoom = 1
    end
    if newfig == true
        fig = figure(title,figsize=(10,10))
        clf()
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

function imview2(x1, x2; figtitle="", caption1="", caption2 = "", zoom=1, color="gray")
    if zoom < 1
        zoom = 1
    end
    fig = figure(figtitle,figsize=(10,5))
    clf()
    nx = size(x1,1)
    ny = size(x1,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))

    subplot(1,2,1);
    imshow(x1[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y], cmap=ColorMap(color), interpolation="none");
    title(caption1)
    subplot(1,2,2);
    imshow(x2[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y], cmap=ColorMap(color), interpolation="none");
    title(caption2)
    tight_layout()
end


function imsurf(z; figtitle="", zoom=1, color="coolwarm")
    if zoom < 1
        zoom = 1
    end
    fig = figure(figtitle,figsize=(10,10))
    clf()
    axis("equal")
    nx = size(z,1)
    ny = size(z,2)
    center_x = div(nx+1,2)
    center_y = div(ny+1,2)
    zoom_span_x = max(1, Int(floor(nx/(2*zoom))))
    zoom_span_y = max(1, Int(floor(ny/(2*zoom))))
    zoomed = z[center_x-zoom_span_x+1:center_x+zoom_span_x, center_y-zoom_span_y+1:center_y+zoom_span_y]
    zoomed[isnan.(zoomed)].=0
    x = collect(1:size(zoomed,1))
    y = collect(1:size(zoomed,2))
    surf(x,y,zoomed, facecolors=get_cmap(color)(zoomed/maximum(zoomed)), linewidth=0.25, rstride=4, cstride=4,antialiased=false, shade=false)
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

function imview4(x1,x2,x3,x4;color="gray")
im = figure("Reconstruction",figsize=(16,4))
clf();
subplot(1,4,1);
imshow(x1, interpolation="none",  cmap=ColorMap(color))
xlabel("Object");
subplot(1,4,2);
imshow(x2, interpolation="none",  cmap=ColorMap(color))
xlabel("Convolved object");
subplot(1,4,3);
imshow(x3, interpolation="none",  cmap=ColorMap(color))
xlabel("Data");
subplot(1,4,4);
imshow(x4, interpolation="none",  cmap=ColorMap(color))
xlabel("Reconstructed");
tight_layout();
end

function imview3(xtruth,xnoisy,xreconst;figtitle="", color="gray")
    if(ndims(xtruth)==1)
        nx = Int(sqrt(length(xtruth)));
        xtruth = reshape(xtruth,nx,nx);
        xnoisy = reshape(xnoisy,nx,nx);
        xreconst = reshape(xreconst,nx,nx);
    end

    fig = figure(figtitle,figsize=(15,5))
    subplot(1,3,1)
    imshow(xtruth, cmap=ColorMap(color), interpolation="none");
    title("Truth")
    subplot(1,3,2)
    imshow(xnoisy, cmap=ColorMap(color), interpolation="none");
    title("Noisy")
    subplot(1,3,3)
    imshow(xreconst, cmap=ColorMap(color), interpolation="none");
    title("Reconstructed")
    tight_layout()
end
