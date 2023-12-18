using Plots
using Measures

"""
    plot_graph_3d(Γ::Union{EquilateralMetricGraph,MetricGraph}; save_as=false, set_title=false, color="gray", size=(300,300))

Plot metric graph 'Γ' on 3d-grid.

### Optional Arguments 
- save_as=false: path to save plot
- set_title=false: optional title on plot
- color="gray": color to plot graph

"""
function plot_graph_3d(Γ::Union{EquilateralMetricGraph,MetricGraph}; save_as=false, set_title=false, color="gray", size=(300,300))
    gr();
    plt = plot3d(layout=1, grid=false,size=size, dpi=300, zlim=(-0.1,0.3), ticklabelsize=10)
    ## plot edges
    for e in edges(Γ.G)
        orig = src(e); term = dst(e)
        xe = LinRange(Γ.coords[orig][1],Γ.coords[term][1],100)
        ye = LinRange(Γ.coords[orig][2],Γ.coords[term][2],100)
        ze = zeros(length(xe))
        plot3d!(plt,xe,ye,ze,lw=3,leg=false,color="gray",linewidth=10)
    end
    ## plot vertices
    for v in vertices(Γ.G)
        scatter!(plt,[Γ.coords[v][1]],[Γ.coords[v][2]],[0],color=color,markerstrokecolor=color,markersize=8)
    end

    if set_title != false
        plot3d!(plt,title=set_title)
    end

    if save_as != false
        savefig(save_as)
    end

    return plt 
end

"""
    plot_function_3d(Γ::Union{EquilateralMetricGraph,MetricGraph}, u::Vector{Function}; save_as=false, set_title=false, color_graph="gray", color_func="cornflowerblue", size=(300,300), lw=3)

Plot function 'u' on metric graph 'Γ' on 3d-grid.

### Optional Arguments 
- save_as=false: path to save plot
- set_title=false: optional title on plot
- color_graph="gray": color to plot graph
- color_func="cornflowerblue": color to plot function

"""
function plot_function_3d(Γ::Union{EquilateralMetricGraph,MetricGraph}, u::Vector{Function}; save_as=false, set_title=false, color_graph="gray", color_func="cornflowerblue", size=(300,300), lw=5, grid_off=false)
    gr();
    plt = plot3d(layout=1,size=size,dpi=300,topmargin=-20mm,xticks=false,yticks=false);
    if grid_off == true
        plt=plot3d!(grid=false,showaxis=false,ticks=false,);
    end
    ### plot edges
    for e in edges(Γ.G)
        orig = src(e); term = dst(e)
        xe = LinRange(Γ.coords[orig][1],Γ.coords[term][1],100)
        ye = LinRange(Γ.coords[orig][2],Γ.coords[term][2],100)
        ze = zeros(length(xe))
        plot3d!(plt,xe,ye,ze,lw=lw-2,leg=false,color="gray")
    end
    ### plot vertices
    for v in vertices(Γ.G)
        scatter!(plt,[Γ.coords[v][1]],[Γ.coords[v][2]],[0],color=color_graph,markerstrokecolor=color_graph,markersize=2)
    end
    ## plot function
    for (j,e) in enumerate(edges(Γ.G))
        orig = src(e); term = dst(e)
        xe = LinRange(Γ.coords[orig][1],Γ.coords[term][1],100)
        ye = LinRange(Γ.coords[orig][2],Γ.coords[term][2],100)
        ze = u[j].(LinRange(0,Γ.ℓ,100))
        plot3d!(xe,ye,ze,lw=lw,leg=false,color=color_func,dpi=300)
    end

    if set_title != false
        plot3d!(plt,title=set_title)
    end

    if save_as != false
        savefig(save_as)
    end
    return plt
end

