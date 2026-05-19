"""Constants used throughout the STROOPWAFEL adaptive sampling module."""
ALPHA_IMF = -2.3
SANA_G = -0.55
SANA_ECC = -0.45

# Log-normal natal kick distribution.
# The kick magnitude v [km/s] follows LogNormal(mu, sigma), meaning
# ln(v) ~ Normal(NATAL_KICK_LOG_MU, NATAL_KICK_LOG_SIGMA).
# With mu=5.67, sigma=0.59 the median kick is exp(5.67) ≈ 291 km/s.
NATAL_KICK_LOG_MU    = 5.67   # mean of ln(v_kick / km s⁻¹)
NATAL_KICK_LOG_SIGMA = 0.59   # std dev of ln(v_kick / km s⁻¹)

R_COEFF = [
    [1.71535900,    0.62246212,     -0.92557761,    -1.16996966,    -0.30631491],
    [6.59778800,    -0.42450044,    -12.13339427,   -10.73509484,   -2.51487077],
    [10.08855000,   -7.11727086,    -31.67119479,   -24.24848322,   -5.33608972],
    [1.01249500,    0.32699690,     -0.00923418,    -0.03876858,    -0.00412750],
    [0.07490166,    0.02410413,     0.07233664,     0.03040467,     0.00197741],
    [0.01077422,    0.00000000,     0.00000000,     0.00000000,     0.00000000],
    [3.08223400,    0.94472050,     -2.15200882,    -2.49219496,    -0.63848738],
    [17.84778000,   -7.45345690,    -48.96066856,   -40.05386135,   -9.09331816],
    [0.00022582,    -0.00186899,    0.00388783,     0.00142402,     -0.00007671]
]

R_SOL_TO_AU = 0.00465047
ZSOL = 0.02
MIN_ENTROPY_CHANGE = 0.01

# Minimum value of (1 - rejection_rate) used in every oversampling and
# normalisation calculation.  Caps the oversampling multiplier at ×200
# (= 2 / 0.01) and prevents near-zero denominators from producing
# astronomical array sizes or overflowing float64 normalisation constants.
MIN_ACTIVE_FRACTION = 0.01
