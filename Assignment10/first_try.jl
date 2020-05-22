using Statistics
################################################################################
# Defining Functions
function get_x(e)
    x_t = Array{Float64}(undef, M+T, S)

    x_t[1,:] = e[1,:]
    for i = 2:(M+T)
        x_t[i,:] = 0.5*x_t[i-1,:] + e[i,:] + e[i-1,:]
    end
    return x_t
end

function get_ψ(d)
    ψ = Array{Float64}(undef, M+T)
    ψ[1] = 1

    for j=2:M+T
        ψ[j] = (j-1-d)/j * ψ[j-1]
    end
    return ψ
end

function get_data(x, ψ)

    y_final = Array{Float64}(undef, M+T, S)

    # Generating matrix to sum over
    xmat = Array{Float64}(undef, M+T, size(x)[2])
    ymat = similar(xmat)

    y = Array{Float64}(undef, M+T)
    x_used = Array{Float64}(undef, M+T)
    x_tmp = zeros(M+T)

    for s=1:S
        # Running sample s
        x_used = x[:,s]

        for j=0:S-1
            x_tmp = zeros(M+T)
            x_tmp[j+1:end] = x_used[1:end-j]
            xmat[:,j+1] = x_tmp
        end

        ymat = xmat.*ψ

        for col=1:size(ymat)[1]
            y[col] = sum(ymat[col,:])
        end
        y_final[:,s] = y
        println(s)
    end
    return y_final
end

function Iy(lambdaj, x)
    bigT = length(x)
    prefactor = (1/2*π*bigT)

    firstpart = Array{Float64}(undef, bigT-1)
    secondpart = similar(firstpart)

    for t=1:bigT-1
        firstpart[t] = x[t] * cos(lambdaj*(t-1))
        secondpart[t] = x[t] * sin(lambdaj*(t-1))
    end
    value = (sum(firstpart))^2 + (sum(secondpart))^2
    return value
end

function get_lambdas(α, tslength)
    little_m = floor(Int,tslength^α)
    lambdas = Array{Float64}(undef, little_m)

    for m=1:little_m
        lambdas[m] = (2*π/tslength) * m
    end
    return lambdas, little_m
end

function objf(α, x, d)
    lambda_list, little_m = get_lambdas(α, length(x))

    first_sum = Array{Float64}(undef, little_m)
    second_sum = similar(first_sum)

    for (m, lambda) in enumerate(lambda_list)
        first_sum[m] = Iy(lambda, x) / lambda^(-2*d)
        second_sum[m] = log(lambda)
    end

    firstpart = (1/little_m) * sum(first_sum)
    secondpart = (2*d/little_m) * sum(second_sum)

    value = log(firstpart) - secondpart
    return value
end

function optimize_overd(α, x, ds)
    values = Array{Float64}(undef, length(ds))

    for (i, d) in enumerate(ds)
        values[i] = objf(α, x, d)
    end
    min_d = d_list[argmin(values)]
    return min_d
end

function optimize_overα(αs, x, ds)
    values = Array{Float64}(undef, size(x)[2], length(αs))

    for (i, α) in enumerate(αs)
        for sample=1:size(x)[2]
            values[sample, i] = optimize_overd(α, x[:,sample], ds)
        end
        println(α)
    end
    return values
end
###################################################################
const M=2000
const T=1000
const S=100
const d=0.25
const d_list = LinRange(-0.45, 0.45, 19)
const α_list = LinRange(0.2, 0.8, 7)


e_t = randn(Float64, M+T+1, S)
e_t[1,:] = zeros(S)

x_one = get_x(e_t)
ψ_one = get_ψ(d)

y = get_data(x_one, ψ_one)
y_use = y[M+1:end,:]

df = optimize_overα(α_list, y_use, d_list)
