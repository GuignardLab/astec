import sys

import astec.utils.common as common
import astec.utils.neighborhood as neighborhood
import astec.utils.neighborhood_figure as neighborhoodf

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

########################################################################################
#
#
#
########################################################################################

def neighborhood_process(experiment, parameters):
    proc = "neighborhood_process"
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, neighborhood.NeighborhoodParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    #
    # should we clean reference here?
    #
    neighborhoods = neighborhood.build_neighborhoods(parameters.atlasFiles , parameters,
                                                     time_digits_for_cell_id=time_digits_for_cell_id)

    if True:
        # neighborhoodf.figures_neighborhood_consistency(neighborhoods, parameters)
        neighborhoodf.figures_neighborhood_score(neighborhoods, parameters)

    if parameters.neighborhood_improvement:
       neighborhood.global_score_improvement(neighborhoods, parameters)
