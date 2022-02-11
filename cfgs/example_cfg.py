import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V10_30Dec21',
    qcd_label = 'V10_30Dec21',
    signal_labels = ['V10_30Dec21_m1', 'V10_30Dec21_m1p5', 'V10_30Dec21_m2', 'V10_30Dec21_m3', 'V10_30Dec21_m4p5'],
    sample_type = 'flat',
    tree_name = 'signal_tree',
    qcd_white_list = '20to300',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = '3cat_0_1_5_significance',
    # ----------------------- #

    # ---    SELECTION    --- #
    selection_label = 'baseline_30Dec21',
    # ----------------------- #

    # ---     WEIGHTS     --- #
    add_weight_hlt = False,
    add_weight_pu = True,
    branch_weight_hlt = 'weight_hlt_HLT_Mu9_IP6_A1_6',
    branch_weight_pu = 'weight_pu_qcd_D',
    # ----------------------- #


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###   Plotting Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    # quantities 
    quantities_label = ['all', 'muonId_study_displacedmuon', 'trackId'],
    # what to plot
    plot_CR = False,
    plot_SR = True,
    plot_dataSig = False,
    # plot style
    do_shape = False,
    do_luminorm = True,
    do_stack = False,
    do_log = False,
    plot_ratio = False,
    add_overflow = True, 
    add_CMSlabel = True,
    add_lumilabel = True,
    CMStag = 'Preliminary',


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    Datacard info    ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    # choose analysis strategy
    do_counting = False,
    do_shape_analysis = True,
    do_shape_TH1 = False,
    lumi_target = 41.6,

    # if counting, or shape_TH1, choose bkg estimation method
    do_ABCD = False,
    do_ABCDHybrid = True,
    do_TF = False,
    do_realData = False,
    ABCD_label = 'cos2d_svprob_0p996',
    sigma_mult_window = 2,

    # if shape analysis
    signal_model_label = 'voigtian',
    background_model_label = 'chebychev',
    do_binned_fit = True,
    do_blind = True,
    nbins = 50, 
    plot_pulls = True,

    # further options
    do_categories = True,
    add_Bc = False,
    plot_prefit = True,
    sigma_B = 472.8e9,


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###     Limits info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    run_blind = True,
)







