using Statistics
################################################################################
# Defining Functions
function get_x(e)
    x_t = Array{Float64}(undef, M+T, S)

    x_t[1,:] = e[1,:]
    for i = 2:(M+T)
        x_t[i,:] = 0.5*x_t[i-1,:] + e[i,:] + 0.3*e[i-1,:]
    end
    return x_t
end

function get_ψ(d)
    ψ = Array{Float64}(undef, M+T)
    ψ[1] = 1

    for j=1:M+T-1
        ψ[j+1] = ((j-1-d)/j) * ψ[j]
    end
    return ψ
end

function get_data(x, ψ)

    y_final = Array{Float64}(undef, M+T, S)

    # Generating matrix to sum over
    xmat = Array{Float64}(undef, M+T, M+T)
    ymat = similar(xmat)

    y = Array{Float64}(undef, M+T)
    x_used = Array{Float64}(undef, M+T)
    x_tmp = zeros(M+T)

    for s=1:S
        # Running sample s
        x_used = x[:,s]

        for j=0:M+T-1
            x_tmp = zeros(M+T)
            x_tmp[j+1:end] = x_used[1:end-j]
            xmat[j+1,:] = x_tmp
        end

        ymat = xmat.*ψ

        for col=1:size(ymat)[2]
            y[col] = sum(ymat[:,col])
        end
        y_final[:,s] = y
        println(s)
    end
    return y_final
end

function Iy(λj, x)
    bigT = length(x)
    prefactor = (1/2*π*bigT)

    firstpart = Array{Float64}(undef, bigT-1)
    secondpart = similar(firstpart)

    for t=1:bigT-1
        firstpart[t] = x[t] * cos(λj*(t-1))
        secondpart[t] = x[t] * sin(λj*(t-1))
    end
    value = (sum(firstpart))^2 + (sum(secondpart))^2
    return value
end

function get_λs(α, tslength)
    little_m = floor(Int,tslength^α)
    λs = Array{Float64}(undef, little_m)

    for m=1:little_m
        λs[m] = (2*π/tslength) * m
    end
    return λs, little_m
end

function objf(α, x, d)
    λ_list, little_m = get_λs(α, length(x))

    first_sum = Array{Float64}(undef, little_m)
    second_sum = similar(first_sum)

    for (m, λ) in enumerate(λ_list)
        first_sum[m] = Iy(λ, x) / λ^(-2*d)
        second_sum[m] = log(λ)
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
    values_LW = Array{Float64}(undef, size(x)[2], length(αs))
    values_PR = similar(values_LW)

    for (i, α) in enumerate(αs)
        for sample=1:size(x)[2]
            values_LW[sample, i] = optimize_overd(α, x[:,sample], ds)
            values_PR[sample, i] = PR(α, x[:,sample])
        end
        println(α)
    end
    return values_LW, values_PR
end

function Ij(λj, x)
    value = Iy(λj, x)
    return log(value)
end

function Rj(λj)
    r = 4 * sin((λj/2))^2
    return -log(r)
end

function OLSestimator(x, y)
    res = inv(x'*x)*(x'*y)
    return res
end

function PR(α, x)
    λ_list, little_m = get_λs(α, length(x))

    x_values = Array{Float64}(undef, little_m, 2)
    y_values = Array{Float64}(undef, little_m, 1)

    x_values[:,1] = ones(little_m)

    for (m, λ) in enumerate(λ_list)
        x_values[m,2] = Rj(λ)
        y_values[m] = Ij(λ, x)
    end

    estimtates = OLSestimator(x_values, y_values)
    return estimtates[2]
end

function get_cis(n, α, type)
    """
    Type 1 = LW, Type 2 = PR
    """

    m = floor(T^α)
    z = 1.96
    if type == 1
        sigma = sqrt(m/4)
    else
        sigma = sqrt((m*π^2/24))
    end

    ci_min = d - (z*sigma)/sqrt(n)
    ci_max = d + z*sigma/sqrt(n)

    return (ci_min, ci_max)
end

function calc_share(x, α, type)
    """
    Type 1 = LW, Type 2 = PR
    """
    n = length(x)
    ci_low, ci_high = get_cis(n, α, type)

    positives = sum([(element > ci_low) & (element < ci_high) for element in x])
    value = positives/n
    return value
end

function calc_bias(x)
    bias = d - mean(x)
    return bias
end

function calc_sigma2(x)
    σ2 = var(x)
    return σ2
end

function cov_probs(LW_values, PR_values, αs)
    probs = Array{Float64}(undef, 2, size(LW_values)[2])
    biases = similar(probs)
    vars = similar(probs)

    for (index, df) in enumerate([LW_values, PR_values])
        for i=1:size(LW_values)[2]
            α = αs[i]
            probs[index,i] = calc_share(df[:,i], α, index)
            biases[index,i] = calc_bias(df[:,i])
            vars[index,i] = calc_sigma2(df[:,i])
        end
    end

    return probs, biases, vars
end
################################################################################
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


# Optimizing
df_LW, df_PR = optimize_overα(α_list, y_use, d_list)
# Calculating final results
co_prob, b, sigma = cov_probs(df_LW, df_PR, α_list)
