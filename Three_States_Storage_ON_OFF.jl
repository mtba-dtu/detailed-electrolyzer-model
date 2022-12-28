#************************************************************************
# Power-to-X Modelling and Optimization
# Three States Model with Storage
#************************************************************************
# Packages
using JuMP
using Gurobi
using DataFrames
using XLSX
#************************************************************************
# Model
include("Input.jl")
m3s = Model(Gurobi.Optimizer)
#************************************************************************
# Variables
@variables m3s begin
    m[T] >= 0 # electrictiy market
    e[T,S] >= 0 # electrolyzer consumption for each segment
    e_tot[T] >= 0 # electrolyzer consumption for each segment
    h[T] >= 0 # hydrogen production
    h_d[T] >= 0 # hydrogen production directly to demand
    d[T] >= 0 # hydrogen sold
    c[T] >= 0 # compressor consumption
    s_in[T] >= 0 # hydrogen stored
    s_out[T] >= 0 # hydrogen used from storage
    soc[T] >= 0 # state of charge of storage (kg)
    z_h[T,S], Bin # specific hydrogen production 
    z_on[T], Bin # on electrolyzer
    z_off[T], Bin # off electrolyzer
    z_sb[T], Bin # standby electrolyzer
    z_start[T], Bin # start electrolyzer
end
#************************************************************************
# Objective function
# Max profit
@objective(m3s, Max, 
    sum(m[t] * lambda_M[t,1] 
        + d[t] * lambda_H 
        - z_start[t] * lambda_start for t=T))
#************************************************************************
# Constraints
# Electricity offer
@constraint(m3s, [t=T],
    m[t] == P_W[t] - e_tot[t] - c[t])
# Total electricity consumption
@constraint(m3s, [t=T],
    e_tot[t] == sum(e[t,s] for s=S))
# Hydrogen production
@constraint(m3s, [t=T],
    h[t] == sum(a[s] * e[t,s] + b[s] * z_h[t,s] for s=S))
# Segment min power boundary
@constraint(m3s, [t=T,s=S],
    e[t,s] >= P_segments[segments][s] * C_E * z_h[t,s])
# Segment max power boundary
@constraint(m3s, [t=T,s=1:length(S)-1],
    e[t,s] <= P_segments[segments][s+1] * C_E * z_h[t,s])
# Hydrogen balance
@constraint(m3s, [t=T],
    h[t] == h_d[t] + s_in[t])
# Demand balance 
@constraint(m3s, [t=T],
    d[t] == h_d[t] + s_out[t])  
# Maximum storage output
@constraint(m3s, [t=T],
    s_out[t] <= C_E * eta_full_load)
# Hydrogen min demand (DAILY)
@constraint(m3s, [tt=1:length(TT)],
    sum(d[t] for t=TT[tt]) >= C_D)
# Maximum electrolyzer power
@constraint(m3s, [t=T],
    e_tot[t] <= C_E * z_on[t])
# Minimum electrolyzer power
@constraint(m3s, [t=T],
    e_tot[t] >= P_min * z_on[t])
# only one efficiency if on or standby
@constraint(m3s, [t=T],
    z_on[t] == sum(z_h[t,s] for s=S))
# States
@constraint(m3s, [t=T],
    1 == z_on[t] + z_off[t])
# Startup cost
@constraint(m3s, [t=T],
    z_start[t] >= (t > 1 ? z_on[t] - z_on[t-1] : 0))
# Compressor consumption
@constraint(m3s, [t=T],
    c[t] == s_in[t] * P_C)
# Max storage fill
@constraint(m3s, [t=T],
    soc[t] <= C_S)
# SOC
@constraint(m3s, [t=T],
    soc[t] == (t > 1 ? soc[t-1] : 0) - s_out[t] + s_in[t])
#************************************************************************
# Solve
optimize!(m3s)
# Report results
println("-------------------------------------");
if termination_status(m3s) == MOI.OPTIMAL
    println(objective_value(m3s))
    output = DataFrame()
    output[!,:CF] = CF
    output[!,:Wind] = P_W
    output[!,:P] = lambda_M
    output[!,:H2_tresh] = lambda_M / lambda_H
    output[!,:Sold_El] = value.(m)
    #e_list = []
    #for t=T push!(e_list,sum(value.(e[t,s] for s=S))) end
    #output[!,:Elec] = e_list
    output[!,:Elec] = value.(e_tot)
    output[!,:ON] = value.(z_on)
    output[!,:S_UP] = value.(z_start)    
    output[!,:H2_prod] = value.(h)
    eff_list = []
    for t=T push!(eff_list,value.(h[t])/sum(value.(e[t,s] for s=S))) end
    output[!,:Elec_Eff] = eff_list
    output[!,:Comp] = value.(c)
    output[!,:H2_store_in] = value.(s_in)
    output[!,:H2_store_out] = value.(s_out)
    output[!,:SOC] = value.(soc)                
    output[!,:H2_sold] = value.(d)            
    show(output) 
    println("\n\nRESULTS:")
    println("Objective value    = $(round(objective_value(m3s), digits=digs)) EUR")
    println("Electricity        = $(round(sum(value.(m[t]) * lambda_M[t] for t=T), digits=digs)) EUR")
    println("Hydrogen           = $(round(sum(value(d[t]) * lambda_H for t=T), digits=digs)) EUR")
    println("Max SOC            = $(maximum(value.(soc)))")
    else
    println(" No solution")
end
println("\n--------------------------------------");
#************************************************************************
# Output to EXCEL
XLSX.writetable("m3s_ON_OFF_$(segments)_segments.xlsx", output, overwrite=true)