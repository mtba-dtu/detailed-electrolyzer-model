# Power-to-X Modelling and Optimization
#************************************************************************
# Input Data
#************************************************************************
using DataFrames
using XLSX
using PyCall
using Plots
digs = 2 # rounding digits
#************************************************************************
Plot = false # output of plots true or false
#************************************************************************
# Segments
segments = 1
# Input
Scenario = DataFrame(XLSX.readtable("scenario_data_2019.xlsx", "scenario"))
T = collect(1:nrow(Scenario))
# Wind farm
C_W = 104.5# nominal power wind in MW
CF = Scenario[!,:cf]
P_W = CF * C_W # wind production 
lambda_M = Scenario[!,:p] # electricity day ahead market price
# Electrolyzer
C_E = C_W/2 # size in MW
sb_load = 0.01 # standby load = 1%
P_sb = C_E * sb_load # standby load = 1%
min_load = 0.15 # minimum load = 15%
P_min = C_E * min_load # minimum load = 15%
p_cell = 30 # cell pressure bar
T_cell = 90 # cell temperature in celsius
i_max = 5000 # maximum cell current density in A/m2
A_cell = 0.2 # cell area in m2
start_cost = 50 # starting cost of production = 50 EUR/MW
lambda_start = C_E * start_cost # starting cost of production
eta_full_load = 17.547 # constant production efficiency kg/MWh
TSO_tariff = 15.06 # TSO grid tariff in EUR/MWh 
lambda_M_in = Scenario[!,:p] .+ TSO_tariff # electricity day ahead market price
# Hydrogen market
lambda_H = 2.10 # EUR per kg
# Hydrogen storage
C_S = C_E * eta_full_load * 24 # max size in kg
soc_0 = 0 # initial storage in MWh
# Compressor
eta_C = 0.75 # mechanical efficiency in %
p_in = 30 # inlet pressure in bar
p_out = 200 # outlet pressure in bar
gamma = 1.4 # adiabatic exponent
T_in = 40 + 273.15 # inlet temperature in K
R = 8.314 # universal gas constant in J/mol*K
M_H2_kg = 2.0159E-03 # molar mass of H2 in kg/mol
P_C = R * T_in / M_H2_kg * gamma / (gamma-1) * 1 / eta_C * ((p_out/p_in)^((gamma-1)/gamma)-1) * 1E-06 / 3600 # compressor consumption in MWh/kg H2

# Daily hours
TT_daily = []
for x in 1:365
    hour = []
    for y in 1:24
        push!(hour, y+(x-1)*24)
    end
    push!(TT_daily,hour)
end

# Monthly hours
TT_monthly = []
push!(TT_monthly, collect(1:31*24)) # Jan
push!(TT_monthly, collect(31*24+1:59*24)) # Feb
push!(TT_monthly, collect(59*24+1:90*24)) # Mar
push!(TT_monthly, collect(90*24+1:120*24)) # Apr
push!(TT_monthly, collect(120*24+1:151*24)) # Mai
push!(TT_monthly, collect(151*24+1:181*24)) # Jun
push!(TT_monthly, collect(181*24+1:212*24)) # Jul
push!(TT_monthly, collect(212*24+1:243*24)) # Aug
push!(TT_monthly, collect(243*24+1:273*24)) # Sep
push!(TT_monthly, collect(273*24+1:304*24)) # Oct
push!(TT_monthly, collect(304*24+1:334*24)) # Nov
push!(TT_monthly, collect(334*24+1:365*24)) # Dec

# Demand
C_D_daily = C_E * eta_full_load * 4 # in kg per day
C_D_monthly = C_D_daily * 30 # in kg per year
C_D_annual = C_D_daily * 365 # in kg per year

# Set demand style
C_D = C_D_daily
TT = TT_daily
#************************************************************************
# Coefficients
a1	= 1.5184
a2	= 1.5421E-03
a3	= 9.523E-05
a4	= 9.84E-08
r1 = 4.45153E-05
r2 = 6.88874E-09
d1 = -3.12996E-06
d2 = 4.47137E-07
s = 0.33824
t1 = -0.01539
t2 = 2.00181
t3 = 15.24178
B1 = 4.50E-05
B2 = 1.02116
B3 = -247.26
B4 = 2.06972
B5 = -0.03571
f11 = 478645.74
f12 = -2953.15
f21 = 1.0396
f22 = -0.00104
F_const = 96485.3321
M_H2 = 2.0159 # molar mass of H2 in kg/mol
HHV = 39.41
#************************************************************************
#Functions
# Reversible cell voltage
function U_rev(Temp)
    Temp_K = Temp + 273.15
    U_rev = a1 - a2 * Temp_K + a3 * Temp_K * log(Temp_K) + a4 * Temp_K^2
    return U_rev
end
# Real cell voltage
function U_cell(Temp,p,i)
    U_cell = U_rev(Temp) + ((r1 + d1) + r2 * Temp + d2 * p) * i + s * log(10,(t1 + t2 / Temp + t3 / Temp^2) * i + 1) 
    return U_cell
end
# Cell power consumption
function P_cell(Temp,p,i)
    P_cell = i * U_cell(Temp,p,i)
    return P_cell
end
# Faraday efficiency (5-parameter)
function eta_F_5(Temp,i)
    eta_F = B1 + B2 * exp((B3 + B4 * Temp + B5 * Temp^2) / i)
    return eta_F
end
# Faraday efficiency
function eta_F(Temp,i)
    eta_F = (i^2 / (f11 + f12 * Temp + i^2)) * (f21 + f22 * Temp)
    return eta_F
end
# Cell production
function M_H_cell(Temp,i)
    M_H_cell = (eta_F(Temp,i) * M_H2 * i) / (2 * F_const)
    M_H_cell_kg_h = M_H_cell * 3.6
    return M_H_cell_kg_h
end
# System production
function M_H_sys(Temp,i,I,n_c)
    M_H_cell = (eta_F(Temp,i) * n_c * M_H2 * I) / (2 * F_const)
    M_H_cell_kg_h = M_H_cell * 3.6
    return M_H_cell_kg_h
end
# Cell efficiency
function eta_cell(Temp,p,i)
    eta_cell = M_H_cell(Temp,i) * HHV / P_cell(Temp,p,i)
    return eta_cell
end
# Number of cells
function n_cell(i_max,A_cell,C_E,Temp,p)
    I_max_cell = i_max * A_cell
    U_max_cell = U_cell(Temp,p,i_max)
    P_max_cell = I_max_cell * U_max_cell
    #n_cell = ceil((C_E * 1000000) / P_max_cell)
    n_cell = (C_E * 1000000) / P_max_cell
    return n_cell
end
# Production curve
function P_curve(P_list,S,i_max,A_cell,C_E,Temp,p)
    N = collect(1:length(P_list))
    a = []
    b = []
    n_x = []
    n_y = []
    n_c = n_cell(i_max,A_cell,C_E,Temp,p)
    include("find_current.jl")
    i_list = py"find_i_from_p"(P_list,C_E,n_c,A_cell,Temp,p)

    for n=N 
        i = i_list[n]
        I = i_list[n] * A_cell
        U = U_cell(Temp,p,i) * n_c
        P = I * U / 1000000
        M = M_H_sys(Temp,i,I,n_c)
        push!(n_x, P)
        push!(n_y, M)        
    end
    for s=S # a and b for segment
        a_s = (n_y[s] - n_y[s+1]) / (n_x[s] - n_x[s+1])
        b_s = n_y[s] - (a_s * n_x[s])
        push!(a, a_s)
        push!(b, b_s)
    end

    return a,b,n_x,n_y
end
#************************************************************************
# Production curve
P_E_min = min_load # minimum load for production
P_E_opt = 0.28231501
P_E_max = 1
P_segments = [
    [P_E_min,
        P_E_max], #1
    [P_E_min,
        P_E_opt,
        P_E_max], #2
    [],
    [P_E_min,
        (P_E_min+P_E_opt)/2,
        P_E_opt,
        (P_E_opt+P_E_max)/2,
        P_E_max], #4
    [],[],[],
    [P_E_min,
        (P_E_min+(P_E_min+P_E_opt)/2)/2,
        (P_E_min+P_E_opt)/2,
        ((P_E_min+P_E_opt)/2+P_E_opt)/2,
        P_E_opt,
        (P_E_opt+(P_E_opt+P_E_max)/2)/2,
        (P_E_opt+P_E_max)/2,
        ((P_E_opt+P_E_max)/2+P_E_max)/2,
        P_E_max], #8
    [],[],[],
    [P_E_min,
        (P_E_min+(P_E_min+P_E_opt)/2)/2,
        (P_E_min+P_E_opt)/2,
        ((P_E_min+P_E_opt)/2+P_E_opt)/2,
        P_E_opt,
        (P_E_opt+(P_E_opt+(P_E_opt+P_E_max)/2)/2)/2,
        (P_E_opt+(P_E_opt+P_E_max)/2)/2,
        ((P_E_opt+(P_E_opt+P_E_max)/2)/2+(P_E_opt+P_E_max)/2)/2,
        (P_E_opt+P_E_max)/2,
        ((P_E_opt+P_E_max)/2+((P_E_opt+P_E_max)/2+P_E_max)/2)/2,
        ((P_E_opt+P_E_max)/2+P_E_max)/2,
        (((P_E_opt+P_E_max)/2+P_E_max)/2+P_E_max)/2,
        P_E_max] #12
]
P_list = P_segments[segments]
S = collect(1:length(P_list)-1)
curve = P_curve(P_list,S,i_max,A_cell,C_E,T_cell,p_cell)
a = curve[1]
b = curve[2]
#************************************************************************
# Non-linear plot
n_c = n_cell(i_max,A_cell,C_E,T_cell,p_cell)
df = DataFrame(i=1:i_max,u=zeros(),eta_F=zeros(),load=zeros(),power=zeros(),prod=zeros(),eff=zeros())
for i in df.i
    I = i * A_cell
    df[i,:u] = U_cell(T_cell,p_cell,i)
    df[i,:eta_F] = eta_F(T_cell,i)
    df[i,:power] = i * A_cell * df[i,:u] * n_c / 1000000
    df[i,:load] = df[i,:power] / C_E 
    df[i,:prod] = M_H_sys(T_cell,i,I,n_c)
    df[i,:eff] = df[i,:prod] / (df[i,:power])
end
#************************************************************************
# Plotting
if Plot == true
# Plot production curve
plot(df[!,:power],df[!,:prod],
        label="Curve",
        line=4,
        color=:red)
plot!(curve[3],curve[4],
    label="Segments",
    line=3,
    color=:blue) 
display(plot!(
    xlims=(3,20),
    xticks=0:2:20,
    ylims=(minimum(curve[1]),maximum(df[!,:prod])*1.1),
    yticks=0:25:350,
    xlabel="Consumption [MWh]",
    ylabel="Production [kg]",
    title="Production Curve",
    legend=:bottomright,
    framestyle =:box)) 
end
