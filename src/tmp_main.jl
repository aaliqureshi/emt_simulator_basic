# d0 = load_results("results0/ieee39_fault.json")
# d1 = load_results("results_line/ieee39_fault.json")
# d2 = load_results("results_full/ieee39_fault.json")

d0 = load_results("results0/ieee14_fault_barq.json")
d1 = load_results("results_line/ieee14_fault_barq.json")
d2 = load_results("results_full/ieee14_fault_barq.json")

using Plots

plotly()
default(label="")
# d0

x = d0[:time]
plot(d0[:time], d0[:delta][:,:]')
plot(d0[:time], d0[:omega][:,:]')

v0_complex = d0[:balance_d] + 1im * d0[:balance_q]
v0 = real.(v0_complex)

v1_complex = d1[:balance_d] + 1im * d1[:balance_q]
v1 = real.(v1_complex)

v2_complex = d2[:balance_d] + 1im * d2[:balance_q]
v2 = real.(v2_complex)


plot(d0[:time], v0')

plot(d1[:time], v1')

plot(d2[:time], v2')

plot(d0[:time], v0' - v1')
plot(d0[:time], v0' - v2')

plot(d0[:time], v1' - v2')

diff1 = v0' - v1'

plot(x, diff1[:,12])

plot(x, v0[1,:])
plot!(x, v1[1,:])
plot!(x, v2[1,:])


### current
i0_complex = d0[:line_id] + 1im * d0[:line_iq]
i0 = real.(i0_complex)

i1_complex = d1[:line_id] + 1im * d1[:line_iq]
i1 = real.(i1_complex)

i2_complex = d2[:line_id] + 1im * d2[:line_iq]
i2 = real.(i2_complex)

plot(d0[:time], i0')
plot(d1[:time], i1')
plot(d2[:time], i2')


plot(x, i0[1,:])
plot!(x, i1[1,:])
plot!(x, i2[1,:])

diff1 = i0' - i1'
diff2 = i0' - i2'
diff3 = i1' - i2'

plot(x, diff1)
plot(x, diff2)
plot(x, diff3)

gen_ids = collect(1:9)
plot(x, d0[:omega][:,:]', label=["gen $i" for i in gen_ids'])

Z = inv(sys.Y)
z_path = Z[models.fault.bus[1], models.generator.idx]

plot(x, v0[models.generator.bus[models.generator.idx[6]][1],:])
plot!(x, v1[models.generator.bus[models.generator.idx[6]][1],:])


plot(x, v0[models.generator.bus[models.generator.idx[8]][1],:])
plot!(x, v1[models.generator.bus[models.generator.idx[8]][1],:])

d1 = v0[models.generator.bus[models.generator.idx[6]][1], :] - v1[models.generator.bus[models.generator.idx[6]][1], :]
d2 = v0[models.generator.bus[models.generator.idx[8]][1], :] - v1[models.generator.bus[models.generator.idx[8]][1], :]

plot(d1)
plot!(d2)

z_path_idx_perm = sortperm(abs.(real.(z_path)))

d1 = v0[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :] - v1[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :]
d2 = v0[models.generator.bus[models.generator.idx[z_path_idx_perm[end]]][1], :] - v1[models.generator.bus[models.generator.idx[z_path_idx_perm[end]]][1], :]

plot(d1)
plot!(d2)



plot(x, d0[:omega][models.generator.idx[z_path_idx_perm], :]', label=["gen $i" for i in models.generator.idx[z_path_idx_perm]'])
dd = v0[models.generator.bus[models.generator.idx[z_path_idx_perm]], :] - v1[models.generator.bus[models.generator.idx[z_path_idx_perm]], :]
plot(dd[:,:]', label=["gen $i" for i in models.generator.idx[z_path_idx_perm]'])


d1 = v0[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :] - v1[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :]
d2 = v0[models.generator.bus[models.generator.idx[z_path_idx_perm[end-7]]][1], :] - v1[models.generator.bus[models.generator.idx[z_path_idx_perm[end-7]]][1], :]

plot(d1)
plot!(d2)


d1 = i0[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :] - i1[models.generator.bus[models.generator.idx[z_path_idx_perm[1]]][1], :]
d2 = i0[models.generator.bus[models.generator.idx[z_path_idx_perm[end-2]]][1], :] - i1[models.generator.bus[models.generator.idx[z_path_idx_perm[end-2]]][1], :]

plot(d1)
plot!(d2)