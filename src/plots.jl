using PyPlot, Statistics


function imdisp(x; colormap="gist_gray")
    imshow(x, ColorMap(colormap))
end

function plot_materials_avg(object_true, x_sol, λ, template)
    N = size(object_true, 1)
    object_template = Int.(crop_to(readfits(template), (N, N)));
    segments = sort(unique(object_template));
    segindx  = [findall( object_template .== segments[n]) for n=1:length(segments)];
    fig = figure(1,(21,3))
    fig.clear();
    axes = fig.subplots(1,7)
    subplot(1,7,1)
    plot(λ*1e9, vec(mean(object_true[segindx[2],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[2],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,2)
    plot(λ*1e9, vec(mean(object_true[segindx[3],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[3],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,3)
    plot(λ*1e9, vec(mean(object_true[segindx[4],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[4],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,4)
    plot(λ*1e9, vec(mean(object_true[segindx[5],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[5],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,5)
    plot(λ*1e9, vec(mean(object_true[segindx[6],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[6],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,6)
    plot(λ*1e9, vec(mean(object_true[segindx[7],:], dims=1)), label="Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[7],:], dims=1)), label="Recontructed")
    legend()
    subplot(1,7,7)
    plot(λ*1e9, vec(mean(object_true[segindx[1],:], dims=1)), label="Background Truth")
    plot(λ*1e9, vec(mean(x_sol[segindx[1],:], dims=1)), label="Recontructed")
    legend()
    tight_layout()
end

function plot_materials_pix(object_true, x_sol, σ_sol, λ, template)
    N = size(object_true, 1)
    object_template = Int.(crop_to(readfits(template), (N, N)));
    allpix = findall(object_template.>0)
    reorder = sortperm(object_template[allpix])
    allpix = allpix[reorder]
    colors=["b", "c", "g", "y", "r", "k"]
    nplots = length(allpix)
    nrows = ceil(Int, sqrt(nplots))
    fig = figure(1,(nrows*2,nrows*2))
    fig.clear();
    axes = fig.subplots(nrows,nrows)
    tight_layout()
    for i=1:nplots
        subplot(nrows,nrows,i)
        indx = allpix[i]
        plot(λ*1e9, object_true[indx,:], label="Truth", color=colors[object_template[indx]], linestyle="solid")
        plot(λ*1e9, x_sol[indx,:], label="Truth", color=colors[object_template[indx]], linestyle="dotted")
        errorbar(λ*1e9,x_sol[indx,:], yerr=σ_sol[indx,:],color =colors[object_template[indx]], ecolor=colors[object_template[indx]], elinewidth=1.0, linestyle="dotted")
    end
    for i=nplots+1:nrows*nrows
        axes[div(i-1,nrows)+1, mod(i-1, nrows)+1].set_visible(false) 
    end
end

function plot_materials_pix_shape(object_true, x_sol, σ_sol, λ, template; ignore=6, layout="vertical")
    N = size(object_true, 1)
    object_template = Int.(crop_to(readfits(template), (N, N)));
    allpix = findall( (object_template.>0) .& (object_template .!= ignore))
    x=[allpix[i][1] for i=1:length(allpix)]
    y=[allpix[i][2] for i=1:length(allpix)]
    ncols = maximum(x)-minimum(x)+1
    nrows = maximum(y)-minimum(y)+1
    colors=["b", "c", "g", "y", "r", "k"]
    nplots = length(allpix)
    
    PyDict(pyimport("matplotlib")."rcParams")["font.family"]=["serif"]
    PyDict(pyimport("matplotlib")."rcParams")["font.size"]=[6]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.major.size"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.major.size"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.major.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.major.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["lines.markeredgewidth"]=[1]
    
    if layout == "vertical"
    fig = figure(1,(ncols, nrows))
    fig.clear();
    axes = fig.subplots(nrows, ncols, sharex=true)
    for ix=minimum(x):maximum(x)
        for iy=minimum(y):maximum(y)
            if object_template[ix,iy] ==0
                axes[iy-minimum(y)+1,ix-minimum(x)+1].set_visible(false) 
            else
            #i= (iy-minimum(y)+1)+(ix-minimum(x))*nrows
            i= (ix-minimum(x)+1)+(iy-minimum(y))*ncols
            subplot(nrows, ncols, i)
            #plot(rand(20),color=colors[object_template[ix,iy]])
            plot(λ*1e9, object_true[ix,iy,:], label="Truth", color=colors[object_template[ix,iy]], linestyle="solid")
            plot(λ*1e9, x_sol[ix,iy,:], label="Truth", color=colors[object_template[ix,iy]], linestyle="dotted")
            errorbar(λ*1e9,x_sol[ix,iy,:], yerr=σ_sol[ix,iy,:],color =colors[object_template[ix,iy]], ecolor=colors[object_template[ix,iy]], elinewidth=1.0, linestyle="dotted")   
            end
            xticks([500, 700, 900])
        end

    end
   tight_layout(pad=0.2)
     
    elseif layout == "horizontal"
    fig = figure(1,(nrows, ncols))
    fig.clear();
    axes = fig.subplots(ncols, nrows, sharex=true)
    for ix=minimum(x):maximum(x)
        for iy=minimum(y):maximum(y)
            if object_template[ix,iy] ==0
                axes[ix-minimum(x)+1,iy-minimum(y)+1].set_visible(false) 
            else
            i= (iy-minimum(y)+1)+(ix-minimum(x))*nrows
            subplot(ncols, nrows, i)
            #plot(rand(20),color=colors[object_template[ix,iy]])
             plot(λ*1e9, object_true[ix,iy,:], label="Truth", color=colors[object_template[ix,iy]], linestyle="solid")
             plot(λ*1e9, x_sol[ix,iy,:], label="Truth", color=colors[object_template[ix,iy]], linestyle="dotted")
             errorbar(λ*1e9,x_sol[ix,iy,:], yerr=σ_sol[ix,iy,:],color =colors[object_template[ix,iy]], ecolor=colors[object_template[ix,iy]], elinewidth=1.0, linestyle="dotted")   
             xticks([500, 700, 900])
            end
        end
    end
    tight_layout(pad=0.2)
    end
end
