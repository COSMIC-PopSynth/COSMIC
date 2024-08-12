from cosmic import _evolvebin
from .filter import parse_column_filters
import operator
import numpy

CHECKSTATE_COLUMNS = numpy.array(
    [
        "binstate",
        "evol_type",
        "mass_1",
        "mass_2",
        "kstar_1",
        "kstar_2",
        "sep",
        "porb",
        "ecc",
        "rrlo_1",
        "rrlo_2",
        "aj_1",
        "aj_2",
        "tms_1",
        "tms_2",
        "massc_1",
        "massc_2",
        "rad_1",
        "rad_2",
        "mass0_1",
        "mass0_2",
        "lum_1",
        "lum_2",
        "radc_1",
        "radc_2",
        "menv_1",
        "menv_2",
        "renv_1",
        "renv_2",
        "omega_spin_1",
        "omega_spin_2",
        "b0_1",
        "b0_2",
        "bacc_1",
        "bacc_2",
        "tacc_1",
        "tacc_2",
        "epoch_1",
        "epoch_2",
        "bhspin_1",
        "bhspin_2",
        "teff_1",
        "teff_2",
    ]
)

DEFAULT_CONDITIONS = [-10e30, -1, 10e30] * CHECKSTATE_COLUMNS.size
DEFAULT_CONDITIONS = numpy.array([DEFAULT_CONDITIONS] * 15)
DEFAULT_DTP_STATE = -1 * numpy.ones(15)


def set_checkstates(timestep_conditions=[]):
    """A function which will detemine different time resolution for different states

    Parameters:
        timestep_conditions : (list) default empty (no dynamic dtp setting)
            a nested list of the many different time resolutions and conditions
            for which you would like to apply those time resolution, e.g.,

            >>> timestep_conditions = [['20.0<mass_1<25.5', '15.5>mass_2>10.0', 'dtp=1.0'],
            >>> ['kstar_1=14', 'lum_1>10.0', 'dtp=0.01'],
            >>> ['2>=binstate>=1', 'dtp=None']]

            The first two sets of conditons would be done in a AND fashion.
            The last condition would end time resolution for the BCM array and it
            would skip to only printing the final state
    """
    # assume that we are not doing any special dtp setting
    _evolvebin.check_dtp.check_dtp = 0

    # Set the default state for checkstate which is that there are no
    # conditional states at which to set a special dtp
    checkstate_array = getattr(_evolvebin.checkstate_array, "checkstate_array")
    checkstate_array[:, :] = DEFAULT_CONDITIONS

    # Again we assume that no condtions exist to set a special dtp
    dtp_state = getattr(_evolvebin.checkstate_params, "dtp_state")
    dtp_state[:] = DEFAULT_DTP_STATE

    # avoid array overflow
    if len(timestep_conditions) > 15:
        raise ValueError("You can only set up to 15 different timestep_conditions")

    for index, condition in enumerate(timestep_conditions):
        # we are checking for conditions
        _evolvebin.check_dtp.check_dtp = 1
        conditions = parse_column_filters(condition)
        for param in conditions:
            # ensure that the param is in the checkstate_array
            if param[0].lower() != "dtp" and param[0].lower() not in CHECKSTATE_COLUMNS:
                raise ValueError(
                    f"`{param[0]}` is not a valid column for timestep_conditions - "\
                     + "valid columns are listed in `cosmic.checkstate.CHECKSTATE_COLUMNS`"
                )

            # find where in the checkstate_array this param is
            param_index = numpy.argwhere(param[0].lower() == CHECKSTATE_COLUMNS)
            if param[0] == "dtp":
                dtp_state = getattr(_evolvebin.checkstate_params, "dtp_state")
                if param[2] == "None":
                    dtp_state[index] = 13700.0
                else:
                    dtp_state[index] = param[2]
                continue

            if param[1] == operator.eq:
                checkstate_array[index, param_index * 3] = param[2]
                checkstate_array[index, param_index * 3 + 2] = param[2]
                checkstate_array[index, param_index * 3 + 1] = 0
            elif param[1] == operator.gt:
                checkstate_array[index, param_index * 3] = param[2]
                checkstate_array[index, param_index * 3 + 1] = 1
            elif param[1] == operator.ge:
                checkstate_array[index, param_index * 3] = param[2]
                checkstate_array[index, param_index * 3 + 1] = 2
            elif param[1] == operator.lt:
                checkstate_array[index, param_index * 3 + 2] = param[2]
                checkstate_array[index, param_index * 3 + 1] = 3
            elif param[1] == operator.le:
                checkstate_array[index, param_index * 3 + 2] = param[2]
                checkstate_array[index, param_index * 3 + 1] = 4
