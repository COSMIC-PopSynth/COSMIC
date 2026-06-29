"""Constants used throughout the STROOPWAFEL adaptive sampling module."""

MIN_ENTROPY_CHANGE = 0.01

# Minimum value of (1 - rejection_rate) used in every oversampling and
# normalisation calculation.  Caps the oversampling multiplier at ×200
# (= 2 / 0.01) and prevents near-zero denominators from producing
# astronomical array sizes or overflowing float64 normalisation constants.
MIN_ACTIVE_FRACTION = 0.01
