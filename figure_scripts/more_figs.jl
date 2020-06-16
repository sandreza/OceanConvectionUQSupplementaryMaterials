using Plots
include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/scripts/utils.jl")
include(pwd() * "/figure_scripts/utils.jl")
case = cases[1]
# pyplot(size = (200,200))
#case = "compromise"
resolution = resolutions[1]
# pyplot()
# get chain
chain, e1, e2 = get_chain(case, resolution[1])
Cᴿ = 0.3
chain[4,:] *= Cᴿ #mutliply by the critical richardson number, for visualizing pdfs
right_bounds[4] *= Cᴿ
###
# http://docs.juliaplots.org/latest/generated/colorschemes/
pyplot(size = (200,200))
save_joint_pdfs = true
save_marginal_pdfs = false
left_bounds .+= 1*1e-4
p = joint_pdfs(chain, left_bounds, right_bounds, parameter_dictionary; indpairs = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], bins = 100, color = cgrad(:blues, rev = false))
#:grays, :inferno, cgrad(:grays, rev = true), cgrad(:ocean, rev = true), :tempo
p[6]

if save_joint_pdfs
    savefig(p[1], pwd() * "/figures/cscn.pdf")
    savefig(p[2], pwd() * "/figures/cscd.pdf")
    savefig(p[3], pwd() * "/figures/csch.pdf")
    savefig(p[4], pwd() * "/figures/cdcn.pdf")
    savefig(p[5], pwd() * "/figures/cnch.pdf")
    savefig(p[6], pwd() * "/figures/cdch.pdf")
end

###
pyplot(size = (300,300))
bins = 100
index = 3
Δx = right_bounds[index]
Δy = 3-0
ratio = 0.25 * Δx / Δy
p3 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = true, ylabel = "pdf", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy))

if save_marginal_pdfs
    savefig(p3, pwd() * "/figures/cd.pdf")
end

index = 4
Δx = right_bounds[index]
Δy = 3-0
ratio = 0.25 * Δx / Δy
p4 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = true, ylabel = "pdf", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy))

if save_marginal_pdfs
    savefig(p4, pwd() * "/figures/ch.pdf")
end

index = 1
Δx = right_bounds[index]
Δy = 3-0
ratio = 0.25 * Δx / Δy
p1 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = true, ylabel = "pdf", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0.0000001, Δy))

if save_marginal_pdfs
    savefig(p1, pwd() * "/figures/cs.pdf")
end

index = 2
Δx = right_bounds[index]
Δy = 3-0
ratio = 0.25 * Δx / Δy
p2 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = true, ylabel = "pdf", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0,Δy))

if save_marginal_pdfs
    savefig(p2, pwd() * "/figures/cn.pdf")
end

###
pyplot(size = (300,300))
bins = 100
index = 3
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p3 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = false, ylabel = "frequency", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy))

if save_marginal_pdfs
    savefig(p3, pwd() * "/figures/cd_nn.pdf")
end

index = 4
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p4 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = false, ylabel = "frequency", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0.0001, Δy))

if save_marginal_pdfs
    savefig(p4, pwd() * "/figures/ch_nn.pdf")
end



index = 1
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p1 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = false, ylabel = "frequency", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0.0000001, Δy))

if save_marginal_pdfs
    savefig(p1, pwd() * "/figures/cs_nn.pdf")
end

index = 2
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p2 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = false, ylabel = "frequency", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy))

if save_marginal_pdfs
    savefig(p2, pwd() * "/figures/cn_nn.pdf")
end


###
pyplot(size = (300,300))
bins = 100
index = 3
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p3 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy), ticks = false)

if save_marginal_pdfs
    savefig(p3, pwd() * "/figures/cd_nn.pdf")
end

index = 4
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p4 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0.0001, Δy), ticks = false)

if save_marginal_pdfs
    savefig(p4, pwd() * "/figures/ch_nn.pdf")
end



index = 1
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p1 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]), bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0.0000001, Δy), ticks = false)

if save_marginal_pdfs
    savefig(p1, pwd() * "/figures/cs_nn.pdf")
end

index = 2
Δx = right_bounds[index]
Δy = 10^5-0
ratio = 1/3 * Δx / Δy
p2 = histogram(chain[index, 1:end-1], xlims = (left_bounds[index], right_bounds[index]),  bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, ylims = (0, Δy), ticks = false)

if save_marginal_pdfs
    savefig(p2, pwd() * "/figures/cn_nn.pdf")
end
