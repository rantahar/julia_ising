using Random
using Statistics
using ArgParse

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    "--Temperature", "-T"
        help = "Temperature parameter in units of the Boltzman constant."
        arg_type = Float64
        default = 0.5
    "--lattice_size", "-L"
        help = "The number of sites in the lattice in each dimension"
        arg_type = Int
        default = 64
    "--measurements", "-m"
        help = "The number of update and measure steps"
        arg_type = Int
        default = 1000
    "--thermalization", "-s"
        help = "Number of thermalization updates"
        arg_type = Int
        default = 1000
    "--updates_per_measurement", "-u"
        help = "Number of updates between measurements"
        arg_type = Int
        default = 4096
end

parsed_args = parse_args(ARGS, arg_settings)

Temperature = parsed_args["Temperature"]
L = parsed_args["lattice_size"]
measurements = parsed_args["measurements"]
thermalization_steps = parsed_args["thermalization"]
updates_per_measurement = parsed_args["updates_per_measurement"]

# Initialize the lattice
spins = ones(L,L)
rng = MersenneTwister(1435)

"Update the lattice a number of times"
function update!(spins, N)
    for i in 1:N
        x = floor(Int, rand(rng, 1)[1]*L) + 1
        y = floor(Int, rand(rng, 1)[1]*L) + 1
        energy_plus = -1*( spins[x, mod1(y+1, L)] + spins[x, mod1(y-1, L)] +
                           spins[mod1(x+1, L), y] + spins[mod1(y-1, L), y] )
        energy_minus = 1*( spins[x, mod1(y+1, L)] + spins[x, mod1(y-1, L)] +
                           spins[mod1(x+1, L), y] + spins[mod1(y-1, L), y] )
        P_plus  = exp( -energy_plus/Temperature )
        P_minus = exp( -energy_minus/Temperature )
        probability_plus = P_plus / (P_plus + P_minus)
        if rand(rng, 1)[1] < probability_plus
            spins[x, y] = 1
        else
            spins[x, y] = -1
        end
    end
end

function calc_energy(spins)
    energy = 2*L*L
    for x in 1:L
        for y in 1:L
            energy -= spins[x, y] * spins[x, mod1(y+1, L)]
            energy -= spins[x, y] * spins[mod1(x+1, L), y]
        end
    end
    return energy
end

function calc_magnetization(spins)
    magnetization = 0
    for x in 1:L
        for y in 1:L
            magnetization += spins[x, y]
        end
    end
    return magnetization
end


update!(spins, thermalization_steps*updates_per_measurement)

E = []
M = []
for i in 1:measurements
    update!(spins, updates_per_measurement)
    push!(E, calc_energy(spins))
    push!(M, calc_magnetization(spins))
end

println("Temperature, Lattice size, Energy, Std of Energy, Magnetization, Std of Magnetization ")
println(Temperature, " ", L, " ", mean(E), " ", std(E), " ", mean(M), " ", std(M))
