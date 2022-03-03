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

    ################################################################################
    #
    # figure file building
    #
    # all figures may have to be generated or only some of them
    #
    generate_figure = (isinstance(parameters.generate_figure, bool) and parameters.generate_figure) or \
                      (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'all') or \
                      (isinstance(parameters.generate_figure, list) and 'all' in parameters.generate_figure)

    #
    # cell-to-cell distance between successive cells in a branch with respect to distance from first cell
    # (from first time point/division to last time point/division)
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'cell-distance-along-branch') \
            or (isinstance(parameters.generate_figure, list) and 'cell-distance-along-branch' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate cell distance along branch file", 1)
        ucontactf.figures_distance_along_branch(atlases, parameters,
                                                time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("... done", 1)

    #
    # plot cell number wrt time without and with temporal registration
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'cell-number-wrt-time') \
            or (isinstance(parameters.generate_figure, list) and 'cell-number-wrt-time' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate cell number wrt time file", 1)
        ucontactf.figures_temporal_registration(atlases, parameters, time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("... done", 1)

    #
    # cell neighbors number wrt total number of cells in the embryo
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'neighbors-wrt-cell-number') \
            or (isinstance(parameters.generate_figure, list) and 'neighbors-wrt-cell-number' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate neighbors histogram file", 1)
        ucontactf.figures_neighbor_histogram(parameters.atlasFiles, parameters,
                                             time_digits_for_cell_id=time_digits_for_cell_id)
        monitoring.to_log_and_console("... done", 1)

    #
    # draw histograms of both right pairing and wrong pairing
    # 2D histograms are at division level
    # 1D histograms are at cell (daughter) level
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'distance-histograms') \
            or (isinstance(parameters.generate_figure, list) and 'distance-histograms' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate division distance histogram file", 1)
        ucontactf.figures_distance_histogram(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

    #
    # draw a dendrogram per division where atlases are grouped with similarity between division
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'division-dendrograms') \
            or (isinstance(parameters.generate_figure, list) and 'division-dendrograms' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate division dendrogram figure file", 1)
        ucontactf.figures_division_dendrogram(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

    #
    # plot embryo volume wrt time without and with temporal registration
    #
    if (isinstance(parameters.generate_figure, str) and parameters.generate_figure == 'embryo-volume') \
            or (isinstance(parameters.generate_figure, list) and 'embryo-volume' in parameters.generate_figure) \
            or generate_figure:
        monitoring.to_log_and_console("... generate embryo volume figure file", 1)
        ucontactf.figures_embryo_volume(atlases, parameters)
        monitoring.to_log_and_console("... done", 1)

    # end of figures
    #
    ################################################################################

    #
    # look for daughter that may improve a global score
    # report it in the console/log file
    # as well as in morphonet selection file
    #
    if parameters.division_permutation_proposal:
        ucontacta.division_permutation_proposal(atlases, parameters)

    return atlases
