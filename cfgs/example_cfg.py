import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V09_06Nov21',
    qcd_label = 'V09_06Nov21',
    signal_labels = ['central_V09_06Nov21_m1', 'central_V09_06Nov21_m3', 'central_V09_06Nov21_m4p5'],
    sample_type = 'flat',
    tree_name = 'signal_tree',
    qcd_white_list = '20to300',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = '3cat_0_1_5_significance',
    # ----------------------- #

    # ---    SELECTION    --- #
    selection_label = 'study_Nov21',
    # ----------------------- #

    # ---     WEIGHTS     --- #
    add_weight_hlt = True,
    add_weight_pu = True,
    branch_weight_hlt = 'weight_hlt_HLT_Mu9_IP6_A1_6',
    branch_weight_pu = 'weight_pu_qcd_A',
    # ----------------------- #


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###   Plotting Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    # quantities 
    quantities_label = ['all', 'muonId_study_displacedmuon', 'trackId'],
    # what to plot
    plot_CR = True,
    plot_SR = False,
    plot_dataSig = False,
    # plot style
    do_shape = True,
    do_luminorm = False,
    do_stack = True,
    do_log = False,
    plot_ratio = False,
    add_overflow = False, 
    add_CMSlabel = True,
    CMStag = 'Preliminary',


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    Datacard info    ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ABCD_label = 'cos2d_svprob_0p996',
    do_ABCD = False,
    do_ABCDHybrid = True,
    do_TF = False,
    do_categories = True,
    add_Bc = False,
    lumi_target = 41.6,
    sigma_B = 472.8e9,
    sigma_mult_window = 2,


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###     Limits info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    run_blind = True,
)







