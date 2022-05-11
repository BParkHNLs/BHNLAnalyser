import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V11_24Apr22',
    qcd_label = 'V11_24Apr22',
    signal_labels = ['V11_24Apr22_m1_large', 'V11_24Apr22_m1p5_large', 'V11_24Apr22_m2_large', 'V11_24Apr22_m3_large', 'V11_24Apr22_m4p5_large'],
    sample_type = 'flat',
    tree_name = 'signal_tree',
    qcd_white_list = '20to300',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = 'V11_24Apr22_permass',
    # ----------------------- #

    # ---    SELECTION    --- #
    selection_label = 'baseline_30Dec21',
    # ----------------------- #

    # ---     WEIGHTS     --- #
    add_weight_hlt = True,
    add_weight_pu = True,
    branch_weight_hlt = 'weight_hlt_HLT_Mu9_IP6_A1_6',
    branch_weight_puqcd = 'weight_pu_qcd_D',
    branch_weight_pusig = 'weight_pu_sig_D',
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
    add_overflow = False, 
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
    lumi_target = 41.601,

    # if counting, or shape_TH1, choose bkg estimation method
    do_ABCD = False,
    do_ABCDHybrid = False,
    do_TF = False,
    do_realData = True,
    ABCD_label = 'cos2d_svprob_0p996',
    sigma_mult_window = 2,

    # if shape analysis
    signal_model_label = 'voigtian',
    use_discrete_profiling = True,
    ## if not discrete profiling, choose fixed background model
    background_model_label = 'chebychev',
    do_binned_fit = True,
    do_blind = True,
    mass_window_size = 3,
    fit_window_size = 10,
    nbins = 100, 
    plot_pulls = False,
    plot_prefit = True,

    # further options
    do_categories = True,
    add_Bc = False,
    sigma_B = 472.8e9,


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###     Limits info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    run_blind = True,
)







