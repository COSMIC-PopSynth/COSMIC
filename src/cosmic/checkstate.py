from cosmic import _evolvebin
from .filter import parse_column_filters
import operator
import numpy as np

CHECKSTATE_COLUMNS = np.array([
    "binstate", "evol_type", "mass_1", "mass_2", "kstar_1", "kstar_2", "sep", "porb", "ecc",
     "rrlo_1", "rrlo_2", "aj_1", "aj_2", "tms_1", "tms_2", "massc_1", "massc_2", "rad_1", "rad_2",
     "mass0_1", "mass0_2", "lum_1", "lum_2", "radc_1", "radc_2", "menv_1", "menv_2", "renv_1", "renv_2",
     "omega_spin_1", "omega_spin_2", "b0_1", "b0_2", "bacc_1", "bacc_2", "tacc_1", "tacc_2",
     "epoch_1", "epoch_2", "bhspin_1", "bhspin_2", "teff_1", "teff_2",
])

# maximum number of timestep_conditions (rows in the Fortran arrays)
N_CONDITIONS = 15

# each row of the checkstate_array holds three slots per column: [lower, flag, upper],
# where the flag is a sum of the following bits (see checkstate.f)
LOWER_STRICT = 1    # lower bound uses > rather than >=
UPPER_STRICT = 2    # upper bound uses < rather than <=
EXCLUDE = 4         # value must not equal the one in checkstate_exclude

DEFAULT_CONDITIONS = [-10e30, 0, 10e30] * CHECKSTATE_COLUMNS.size
DEFAULT_CONDITIONS = np.array([DEFAULT_CONDITIONS] * N_CONDITIONS)
DEFAULT_EXCLUDE = np.zeros((N_CONDITIONS, CHECKSTATE_COLUMNS.size))
DEFAULT_DTP_STATE = -1 * np.ones(N_CONDITIONS)

# which bounds each operator sets and which flag bit it adds
OPERATOR_BOUNDS = {
    operator.eq: (("lower", "upper"), 0),
    operator.gt: (("lower",), LOWER_STRICT),
    operator.ge: (("lower",), 0),
    operator.lt: (("upper",), UPPER_STRICT),
    operator.le: (("upper",), 0),
    operator.ne: (("exclude",), EXCLUDE),
}


def _parse_dtp(value, condition):
    """Convert the value of a `dtp` condition to a float, mapping None to 13700 Myr"""
    if isinstance(value, str):
        if value.lower() == "none":
            return 13700.0
        raise ValueError(f"`dtp={value}` in {condition} is not valid - dtp must be a number or None")
    if value < 0:
        raise ValueError(f"`dtp={value}` in {condition} is not valid - dtp must not be negative")
    return value


def set_checkstates(timestep_conditions=[]):
    """Set the time resolution of the bcm array for different states of the evolution

    Parameters
    ----------
    timestep_conditions : `list`, optional
        A nested list of conditions and the time resolution (``dtp``, in Myr) to use when
        they are met, e.g.

        >>> timestep_conditions = [['20.0<mass_1<25.5', '15.5>mass_2>10.0', 'dtp=1.0'],
        ...                        ['kstar_1=14', 'lum_1>10.0', 'kstar_2!=13', 'dtp=0.01'],
        ...                        ['2>=binstate>=1', 'dtp=None']]

        Each sublist contains exactly one ``dtp`` and any number of conditions, which are
        combined in an AND fashion. Supported operators are ``=``, ``!=``, ``<``, ``<=``,
        ``>`` and ``>=``. A column can have both a lower and upper bound, given either as a
        chained comparison or as two separate conditions, plus one ``!=``. A sublist with
        only a ``dtp`` always applies, and if several sublists apply then the last one is
        used. ``dtp=None`` ends the time resolution so that only the final state is printed.
        Default is an empty list (no dynamic dtp setting).

    Raises
    ------
    ValueError
        If more than 15 sublists are given, a column is not in `CHECKSTATE_COLUMNS`, an
        operator is not supported, the same bound is set twice for a column in one sublist,
        or a sublist does not have exactly one valid ``dtp``
    """
    # avoid array overflow
    if len(timestep_conditions) > N_CONDITIONS:
        raise ValueError(f"You can only set up to {N_CONDITIONS} different timestep_conditions")

    # build the arrays first so that nothing is changed if a condition is invalid
    conditions_array = DEFAULT_CONDITIONS.copy()
    exclude_array = DEFAULT_EXCLUDE.copy()
    dtp_array = DEFAULT_DTP_STATE.copy()

    for index, condition in enumerate(timestep_conditions):
        has_dtp = False
        bounds_set = set()
        for name, op, value in parse_column_filters(condition):
            name = name.lower()
            if name == "dtp":
                if op != operator.eq:
                    raise ValueError(f"dtp must be set with `=` (e.g. 'dtp=0.5') in {condition}")
                if has_dtp:
                    raise ValueError(f"dtp is set more than once in {condition}")
                dtp_array[index] = _parse_dtp(value, condition)
                has_dtp = True
                continue

            # ensure that the param is in the checkstate_array
            if name not in CHECKSTATE_COLUMNS:
                raise ValueError(
                    f"`{name}` is not a valid column for timestep_conditions - "
                    + "valid columns are listed in `cosmic.checkstate.CHECKSTATE_COLUMNS`"
                )
            if isinstance(value, str):
                raise ValueError(f"`{value}` is not a valid value for `{name}` in {condition}")
            if op not in OPERATOR_BOUNDS:
                raise ValueError(
                    f"Unsupported operator for `{name}` in {condition} - "
                    + "use one of =, !=, <, <=, >, >="
                )

            # each bound can only be set once per column in a condition
            bounds, bit = OPERATOR_BOUNDS[op]
            for bound in bounds:
                if (name, bound) in bounds_set:
                    raise ValueError(f"`{name}` has more than one {bound} bound in {condition}")
                bounds_set.add((name, bound))

            # find where in the checkstate_array this param is and set its bounds and flag
            col = int(np.flatnonzero(CHECKSTATE_COLUMNS == name)[0])
            if "lower" in bounds:
                conditions_array[index, 3 * col] = value
            if "upper" in bounds:
                conditions_array[index, 3 * col + 2] = value
            if "exclude" in bounds:
                exclude_array[index, col] = value
            conditions_array[index, 3 * col + 1] = int(conditions_array[index, 3 * col + 1]) | bit

        if not has_dtp:
            raise ValueError(f"No dtp was given in {condition} - each condition needs a dtp")

    _evolvebin.checkstate_array.checkstate_array[:, :] = conditions_array
    _evolvebin.checkstate_exclude.checkstate_exclude[:, :] = exclude_array
    _evolvebin.checkstate_params.dtp_state[:] = dtp_array

    # only check the conditions in the Fortran if there are any
    _evolvebin.check_dtp.check_dtp = int(len(timestep_conditions) > 0)
