import sys

import astec.utils.common as common
import astec.utils.atlas_embryo as uatlase
import astec.utils.atlas_division as uatlasd

#
#
#
#
#

monitoring = common.Monitoring()

_instrumented_ = False

########################################################################################
#
# classes
# - computation parameters
#
########################################################################################


def atlas_embryo_process(experiment, parameters):
    proc = "atlas_embryo_process"
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, uatlase.AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    #
    # read atlas/properties files
    #
    atlases = uatlasd.DivisionAtlases(parameters=parameters)

    #
    # copy from files the properties of interest
    # and temporally register the atlases
    #
    atlases.add_atlases(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    #
    # division based part
    # build division atlas
    #
    atlases.build_division_atlases(parameters)

    #
    # figures (if required)
    #
    atlases.generate_figure(parameters, time_digits_for_cell_id=experiment.get_time_digits_for_cell_id())

    return atlases
