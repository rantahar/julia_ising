using Random

beta = 0.5
L = 64
updates = 1000000
spins = ones(L,L)

rng = MersenneTwister(1435)

function update(spins, N)
    for i in 1:N
        x = floor(Int, rand(rng, 1)[1]*L) + 1
        y = floor(Int, rand(rng, 1)[1]*L) + 1
        energy_plus = -1*( spins[x, mod1(y+1, L)] + spins[x, mod1(y-1, L)] +
                           spins[mod1(x+1, L), y] + spins[mod1(y-1, L), y] )
        energy_minus = 1*( spins[x, mod1(y+1, L)] + spins[x, mod1(y-1, L)] +
                           spins[mod1(x+1, L), y] + spins[mod1(y-1, L), y] )
        P_plus  = exp( -beta*energy_plus )
        P_minus = exp( -beta*energy_minus )
        probability_plus = P_plus / (P_plus + P_minus)
        if rand(rng, 1)[1] < probability_plus
            spins[x, y] = 1
        else
            spins[x, y] = -1
        end
    end
end

function calc_energy(spins)
    energy = 0
    for x in 1:L
        for y in 1:L
            energy -= spins[x, y] * spins[x, mod1(y+1, L)]
            energy -= spins[x, y] * spins[mod1(x+1, L), y]
        end
    end
    return energy
end


update(spins, updates)

println("The energy is ", calc_energy(spins))
