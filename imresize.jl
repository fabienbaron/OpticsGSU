if Pkg.installed("Images.jl") == nothing
           println("Please install the Images package using Pkg.add(\"Images\")\n");
end

using Images
using FITSIO
using PyPlot

function imscale(image, scale)
newsize = (Int(round(scale*size(image,1))), Int(round(scale*size(image,2))))
return Images.imresize(image, newsize)
end

function imresize(image, newsize)
       return Images.imresize(image, newsize)
end

function test_resize()
image=read(FITS("saturn.fits")[1])
#image = image / sum(image)
imshow(image)
clf();
imshow(imscale(image,0.5))


end
