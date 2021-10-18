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
    atlases = ucontacta.Atlases()
    atlases.build_neighborhoods(parameters.atlasFiles, parameters, time_digits_for_cell_id=time_digits_for_cell_id)

    if parameters.generate_figure:

        #
        # draw a graph of reference/atlas embryos per division where edge are valued by the sum of scores
        # (do the right pairing between reference/atlas embryos)
        #
        # ucontactf.figures_division_score(atlases.get_neighborhoods(), parameters)

        #
        # draw score histogram of both right pairing and wrong pairing
        # this is then done at the cell level, not at the division level
        #
        # ucontactf.figures_histogram_scores(atlases.get_neighborhoods(), parameters)

        #
        # draw 2D histogram of right/wrong pairing scores of both daugthers
        # this is then done at the division level
        #
        # ucontactf.figures_histogram2D_scores(atlases, parameters)

        #
        # draw a graph per division where edges rae valued with probability
        #
        # ucontactf.figures_division_probability(atlases, parameters)

        #
        # draw a graph per division where edges rae valued with probability
        #
        ucontactf.figures_division_hierarchical_clustering(atlases, parameters)

        if False and parameters.use_common_neighborhood:
            #
            # do a pca analysis of L1-normalized surface contact vectors per cell/neighborhood
            #
            ucontactf.figures_cell_neighborhood_pca(atlases.get_neighborhoods(), parameters)

    if parameters.naming_improvement:
        # ucontacta.global_score_improvement(atlases.get_neighborhoods(), parameters)
        ucontacta.global_probability_improvement(atlases, parameters)
