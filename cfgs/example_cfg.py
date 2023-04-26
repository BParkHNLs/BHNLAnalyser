import sys
sys.path.append('./objects')
from config import Config


cfg = Config(
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###    General Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#

    # ---     SAMPLES     --- #
    data_label = 'V13_06Feb23_fullBPark',
    qcd_label = 'V13_06Feb23',
    signal_labels = [
      'V13_06Feb23_m1',
      'V13_06Feb23_m1p5',
      'V13_06Feb23_m2',
      'V13_06Feb23_m3',
      'V13_06Feb23_m4p5',
      'V13_06Feb23_m5p5',
      'V42_06Feb23_m1p02',
      'V42_06Feb23_m1p04',
      'V42_06Feb23_m1p06',
      'V42_06Feb23_m1p08',
      'V42_06Feb23_m1p1',
      'V42_06Feb23_m1p12',
      #'V42_06Feb23_m1p14',
      'V42_06Feb23_m1p16',
      'V42_06Feb23_m1p18',
      'V42_06Feb23_m1p2',
      'V42_06Feb23_m1p22',
      'V42_06Feb23_m1p24',
      'V42_06Feb23_m1p26',
      'V42_06Feb23_m1p28',
      'V42_06Feb23_m1p3',
      'V42_06Feb23_m1p32',
      'V42_06Feb23_m1p34',
      'V42_06Feb23_m1p36',
      'V42_06Feb23_m1p38',
      'V42_06Feb23_m1p4',
      'V42_06Feb23_m1p42',
      'V42_06Feb23_m1p44',
      'V42_06Feb23_m1p46',
      'V42_06Feb23_m1p48',
      'V42_06Feb23_m1p53',
      'V42_06Feb23_m1p56',
      'V42_06Feb23_m1p59',
      #'V42_06Feb23_m1p62',
      'V42_06Feb23_m1p65',
      #'V42_06Feb23_m1p68',
      'V42_06Feb23_m1p71',
      'V42_06Feb23_m1p74',
      'V42_06Feb23_m1p77',
      'V42_06Feb23_m1p8',
      'V42_06Feb23_m1p83',
      #'V42_06Feb23_m1p86',
      'V42_06Feb23_m1p89',
      'V42_06Feb23_m1p92',
      'V42_06Feb23_m1p95',
      #'V42_06Feb23_m1p98',
      'V42_06Feb23_m2p05',
      'V42_06Feb23_m2p1',
      'V42_06Feb23_m2p15',
      'V42_06Feb23_m2p2',
      'V42_06Feb23_m2p25',
      'V42_06Feb23_m2p3',
      'V42_06Feb23_m2p35',
      'V42_06Feb23_m2p4',
      'V42_06Feb23_m2p45',
      'V42_06Feb23_m2p5',
      'V42_06Feb23_m2p55',
      'V42_06Feb23_m2p6',
      'V42_06Feb23_m2p65',
      'V42_06Feb23_m2p7',
      'V42_06Feb23_m2p75',
      'V42_06Feb23_m2p8',
      'V42_06Feb23_m2p85',
      'V42_06Feb23_m2p9',
      'V42_06Feb23_m2p95',
      'V42_06Feb23_m3p05',
      'V42_06Feb23_m3p1',
      'V42_06Feb23_m3p15',
      'V42_06Feb23_m3p2',
      'V42_06Feb23_m3p25',
      'V42_06Feb23_m3p3',
      'V42_06Feb23_m3p35',
      'V42_06Feb23_m3p4',
      'V42_06Feb23_m3p45',
      'V42_06Feb23_m3p5',
      'V42_06Feb23_m3p55',
      'V42_06Feb23_m3p6',
      'V42_06Feb23_m3p65',
      'V42_06Feb23_m3p7',
      'V42_06Feb23_m3p75',
      'V42_06Feb23_m3p8',
      'V42_06Feb23_m3p85',
      'V42_06Feb23_m3p9',
      'V42_06Feb23_m3p95',
      'V42_06Feb23_m4p0',
      'V42_06Feb23_m4p1',
      'V42_06Feb23_m4p2',
      #'V42_06Feb23_m4p3',
      #'V42_06Feb23_m4p4',
      'V42_06Feb23_m4p6',
      'V42_06Feb23_m4p7',
      'V42_06Feb23_m4p8',
      'V42_06Feb23_m4p9',
      'V42_06Feb23_m5p0',
      'V42_06Feb23_m5p1',
      'V42_06Feb23_m5p2',
      'V42_06Feb23_m5p3',
      'V42_06Feb23_m5p4',
      'V42_06Feb23_m5p5',
      'V42_06Feb23_m5p6',
      'V42_06Feb23_m5p7',
      'V42_06Feb23_m5p8',
      'V42_06Feb23_m5p9',
      'V42_06Feb23_m6p0',
    ],
    sample_type = 'flat',
    tree_name = 'signal_tree',
    qcd_white_list = '20to300',
    # ----------------------- #

    # ---     POINTS      --- #
    ctau_points_label = 'baseline',
    # ----------------------- #

    # ---   REWEIGHTING   --- #
    reweighting_strategy = 'inclusive',
    # ----------------------- #

    # ---    CATEGORIES   --- #
    categories_label = 'categories_0_50_150_Bc',
    # ----------------------- #

    # ---    SELECTION    --- #
    baseline_selection_label = 'baseline_08Aug22',
    do_cutbased = False,
    do_mva = True,

    # if mva selection
    do_parametric = True,
    training_label = 'V13_06Feb23_2023Apr06_14h13m31s',
    cut_score = 0.99,
    # ----------------------- #

    # ---     WEIGHTS     --- #
    add_weight_hlt = True,
    add_weight_pu = True,
    add_weight_muid = True,
    branch_weight_hlt = 'weight_hlt_fullBpark',
    branch_weight_puqcd = 'weight_pu_qcd_D',
    branch_weight_pusig = 'weight_pu_sig_tot',
    branch_weight_mu0id = 'weight_mu0_softid',
    branch_weight_muid = 'weight_mu_looseid',
    # ----------------------- #


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###   Plotting Info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    # quantities 
    quantities_label = ['small'],
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
    lumi_target = -1,

    # if counting, or shape_TH1, choose bkg estimation method
    do_ABCD = False,
    do_ABCDHybrid = False,
    do_TF = False,
    do_realData = True,
    ABCD_label = 'cos2d_svprob_0p996',
    sigma_mult_window = 2,

    # if shape analysis
    signal_model_label = 'doubleCB',
    use_discrete_profiling = True,
    ## if not discrete profiling, choose fixed background model
    background_model_label = 'chebychev',
    do_binned_fit = True,
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
    resolution_p0 = 6.98338e-04, 
    resolution_p1 = 7.78382e-03, 


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ###     Limits info     ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    do_blind = False,


    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    ### Interpretation info ###
    #-#-#-#-#-#-#-#-#-#-#-#-#-#
    muon_eoslabel = '23_04_08_unblinding_Bc_fullscan_nobernstein_v2',
    electron_eoslabel = '07_04_23_fullStat',
    coupling_scenarios_label = 'triangle',
)







