#!/usr/bin/env python3

import os
import time
from argparse import ArgumentParser
import sys

#
# local imports
# add ASTEC subdirectory
#

import astec.utils.common as common
import astec.algorithms.contact_naming as acontactn
import astec.utils.contact_atlas as ucontacta
import astec.utils.properties as properties
import astec.utils.diagnosis as diagnosis
from astec.wrapping.cpp_wrapping import path_to_vt

#
#
#
#
#


def _set_options(my_parser):
    """

    :param my_parser:
    :return:
    """
    proc = "_set_options"
    if not isinstance(my_parser, ArgumentParser):
        print(proc + ": argument is not of type ArgumentParser")
        return
    #
    # common parameters
    #

    my_parser.add_argument('-p', '--parameters',
                           action='store', dest='parameterFile', const=None,
                           help='python file containing parameters definition')
    my_parser.add_argument('-e', '--embryo-rep',
                           action='store', dest='embryo_path', const=None,
                           help='path to the embryo data')

    #
    # other options
    #

    my_parser.add_argument('-write-selection', '--write-selection', '-write-selections', '--write-selections',
                           action='store_const', dest='write_selection',
                           default=False, const=True,
                           help='write out morphonet selection files')

    #
    # control parameters
    #

    my_parser.add_argument('-k', '--keep-temporary-files',
                           action='store_const', dest='keepTemporaryFiles',
                           default=False, const=True,
                           help='keep temporary files')

    my_parser.add_argument('-f', '--force',
                           action='store_const', dest='forceResultsToBeBuilt',
                           default=False, const=True,
                           help='force building of results')

    my_parser.add_argument('-v', '--verbose',
                           action='count', dest='verbose', default=2,
                           help='incrementation of verboseness')
    my_parser.add_argument('-nv', '--no-verbose',
                           action='store_const', dest='verbose', const=0,
                           help='no verbose at all')
    my_parser.add_argument('-d', '--debug',
                           action='count', dest='debug', default=0,
                           help='incrementation of debug level')
    my_parser.add_argument('-nd', '--no-debug',
                           action='store_const', dest='debug', const=0,
                           help='no debug information')

    doc = "print the list of parameters (with explanations) in the console and exit. "
    doc += "If a parameter file is given, it is taken into account"
    my_parser.add_argument('-pp', '--print-param',
                           action='store_const', dest='printParameters',
                           default=False, const=True, help=doc)

    return


#
#
# main function
#
#


def main():

    ############################################################
    #
    # generic part
    #
    ############################################################

    #
    # initialization
    #
    start_time = time.localtime()
    monitoring = common.Monitoring()
    experiment = common.Experiment()

    #
    # reading command line arguments
    # and update from command line arguments
    #
    parser = ArgumentParser(description='Naming')
    _set_options(parser)
    args = parser.parse_args()

    monitoring.update_from_args(args)
    experiment.update_from_args(args)

    if args.printParameters:
        parameters = acontactn.NamingParameters()
        if args.parameterFile is not None and os.path.isfile(args.parameterFile):
            experiment.update_from_parameter_file(args.parameterFile)
            parameters.update_from_parameter_file(args.parameterFile)
        experiment.print_parameters(directories=[])
        parameters.print_parameters()
        sys.exit(0)

    #
    # read input file(s) from args, write output file from args
    #

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    if args.parameterFile is None:
        prop = properties.read_dictionary(args.inputFile, inputpropertiesdict={})
        prop = properties.set_fate_from_names(prop, time_digits_for_cell_id=time_digits_for_cell_id)
        prop = properties.set_color_from_fate(prop)
        if args.outputFile is not None:
            properties.write_dictionary(args.outputFile, prop)
        sys.exit(0)

    #
    # reading parameter files
    # and updating parameters
    #
    parameter_file = common.get_parameter_file(args.parameterFile)
    experiment.update_from_parameter_file(parameter_file)

    #
    # set
    # 1. the working directory
    #    that's where the logfile will be written
    # 2. the log file name
    #    it creates the logfile dir, if necessary
    #
    # experiment.working_dir = experiment.post_dir
    monitoring.set_log_filename(experiment, __file__, start_time)

    #
    # keep history of command line executions
    # and copy parameter file
    #
    experiment.update_history_at_start(__file__, start_time, parameter_file, path_to_vt())
    # experiment.copy_stamped_file(start_time, parameter_file)

    #
    # copy monitoring information into other "files"
    # so the log filename is known
    #
    common.monitoring.copy(monitoring)

    #
    # write generic information into the log file
    #
    # monitoring.write_configuration()
    # experiment.write_configuration()

    # experiment.write_parameters(monitoring.log_filename)

    ############################################################
    #
    # specific part
    #
    ############################################################

    #
    # copy monitoring information into other "files"
    # so the log filename is known
    #
    acontactn.monitoring.copy(monitoring)
    properties.monitoring.copy(monitoring)
    diagnosis.monitoring.copy(monitoring)
    ucontacta.monitoring.copy(monitoring)

    #
    # manage parameters
    # 1. initialize
    # 2. update parameters
    # 3. write parameters into the logfile
    #

    parameters = acontactn.NamingParameters()
    parameters.update_from_parameter_file(parameter_file)
    parameters.write_parameters(monitoring.log_filename)

    #
    # processing
    #

    prop = acontactn.naming_process(experiment, parameters)

    if args.write_selection or parameters.write_selection:
        time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
        properties.write_morphonet_selection(prop, time_digits_for_cell_id=time_digits_for_cell_id,
                                             directory=parameters.outputDir)
    #
    # end of execution
    # write execution time in both log and history file
    #
    end_time = time.localtime()
    monitoring.update_execution_time(start_time, end_time)
    experiment.update_history_execution_time(__file__, start_time, end_time)

    monitoring.to_console('Total computation time = ' + str(time.mktime(end_time) - time.mktime(start_time)) + ' s')


#
#
# main call
#
#


if __name__ == '__main__':
    main()
