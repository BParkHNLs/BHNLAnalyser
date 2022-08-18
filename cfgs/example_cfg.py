import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V12_08Aug22',
    qcd_label = 'V12_08Aug22',
    signal_labels = ['central_V11_24Apr22_benchmark'],
    sample_type = 'flat',
    tree_name = 'signal_tree',
    qcd_white_list = '20to300',
    # ----------------------- #

    # ---     POINTS      --- #
    points_label = 'generated',
    # ----------------------- #

    # ---   REWEIGHTING   --- #
    reweighting_strategy = 'exclusive_fromlargerctau',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = 'inclusive',
    # ----------------------- #

    # ---    SELECTION    --- #
    selection_label = 'baseline_08Aug22',
    # ----------------------- #

    # ---     WEIGHTS     --- #
    add_weight_hlt = True,
    add_weight_pu = True,
    add_weight_muid = True,
    branch_weight_hlt = 'weight_hlt_D1_tag_fired_HLT_Mu9_IP6_or_HLT_Mu12_IP6_ptdxysig_bsrdst_newcat_max3e6_smalltable_v2',
    branch_weight_puqcd = 'weight_pu_qcd_D',
    branch_weight_pusig = 'weight_pu_sig_D',
    branch_weight_mu0id = 'weight_mu0_softid',
    branch_weight_muid = 'weight_mu_looseid',
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
    do_stack = False,
    do_log = False,
    plot_ratio = False,
    add_overflow = False, 
    add_CMSlabel = True,
    add_lumilabel = False,
    CMStag = '"Preliminary"',


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
    plot_prefit = False,

    # further options
    do_categories = True,
    add_Bc = False,
    sigma_B = 472.8e9,
    lhe_efficiency = 0.08244,
    do_tdrstyle = False,


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###     Limits info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    run_blind = True,
)







