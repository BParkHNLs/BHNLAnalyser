import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V07_18Aug21',
    qcd_label = 'V07_18Aug21',
    signal_label = 'private',
    sample_type = 'flat',
    tree_name = 'signal_tree',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = 'standard',
    # ----------------------- #

    # ---    SELECTION    --- #
    selection_label = 'standard',
    # ----------------------- #

    # add weights

    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###   Plotting Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    # quantities
    quantities_label = 'small',
    # what to plot
    plot_CR = True,
    plot_SR = False,
    plot_dataSig = False,
    # plot style
    do_shape = True,
    do_luminorm = False,
    do_stack = False,
    do_log = False,
    add_overflow = False,
    add_CMSlabel = True,
    CMStag = 'Preliminary',


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    Datacard info    ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ABCD_label = 'cos2d_svprob',
    do_categories = True,
    add_Bc = False,

)







