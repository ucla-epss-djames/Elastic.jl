using FFTW

export autocorr_fft, tau_int

# Choose a fast FFT length ≥ n (products of 2,3,5 are usually optimal for FFTW)
nextfast(n::Integer) = Base.nextprod((2,3,5), n)

"""
    autocorr_fft(x; dt=1.0, maxlag=nothing, unbiased=true)

Compute the normalized autocorrelation ρ(τ) of a real-valued time series `x`
via FFT (O(N log N)).

- `dt`      : sampling interval (sets τ units only)
- `maxlag`  : number of lags to return (default: all lags)
- `unbiased`: divide by (N-k) instead of N before normalizing

Returns: (τ::Vector{Float64}, ρ::Vector{Float64})
"""
function autocorr_fft(x; dt=1.0, maxlag=nothing, unbiased=true)
    x = collect(skipmissing(x))         # drop missing values if present
    x = Float64.(x)
    N = length(x)
    N == 0 && error("x is empty")

    # mean-center
    x .-= mean(x)

    # zero-pad to avoid circular wrap-around
    M = nextfast(2N)                    # M ≥ 2N, fast for FFTW
    xpad = zeros(Float64, M)
    @inbounds xpad[1:N] = x

    # Wiener–Khinchin: IFFT of power spectrum = autocovariance
    X = rfft(xpad)
    S = abs2.(X)
    c = irfft(S, M)[1:N]                # linear autocovariance (first N lags)

    # normalization (biased vs unbiased)
    if unbiased
        c ./= (N .- (0:N-1))            # N, N-1, …, 1
    else
        c ./= N
    end

    ρ = c ./ c[1]                       # normalized ACF

    K = maxlag === nothing ? N : min(maxlag, N)
    τ = (0:K-1) .* dt
    return τ, ρ[1:K]
end

# --- quick integral time with a sensible cutoff (stop at first ≤ 0) ---
function tau_int(τ, ρ)
    idx = findfirst(ρ .<= 0)
    K = isnothing(idx) ? length(ρ) : max(idx-1, 1)
    return τ[2] * (1 + 2*sum(ρ[2:K]))   # Δτ*(1+2*sum_{k>=1} ρ_k)
end
