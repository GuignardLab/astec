import sys

import astec.utils.common as common
import astec.utils.contact_atlas as ucontacta
import astec.utils.contact_figure as ucontactf

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


def contact_atlas_process(experiment, parameters):
    proc = "contact_atlas_process"
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ucontacta.AtlasParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    #
    # should we clean reference here?
    #
    atlases = ucontacta.Atlases(parameters=parameters)
    atlases.build_neighborhoods(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    if parameters.generate_figure:

        #
        # draw a graph of reference/atlas embryos per division where edge are valued by the sum of scores
        # (do the right pairing between reference/atlas embryos)
        #
        monitoring.to_log_and_console("... generate division graph figure file", 1)
        ucontactf.figures_division_graph(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

        #
        # draw histograms of both right pairing and wrong pairing
        # 2D histograms are at division level
        # 1D histograms are at cell (daughter) level
        #
        monitoring.to_log_and_console("... generate division distance histogram file", 1)
        ucontactf.figures_distance_histogram(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

        #
        # draw a graph per division where edges are valued with similarity between division
        #
        monitoring.to_log_and_console("... generate division dendrogram figure file", 1)
        ucontactf.figures_division_dendrogram(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

    if parameters.naming_improvement:
        # ucontacta.global_score_improvement(atlases.get_neighborhoods(), parameters)
        ucontacta.division_improvement(atlases, parameters)
