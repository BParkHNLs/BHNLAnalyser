#!/bin/bash

dirlabel='selection_significance_diff_dataD1_v3'
quantity='hnl_cos2d'
action='print_significance'

# list of cuts
# syntax (m1 m3 m4p5)
cut_pi_pt_lxy0to1_OS=("0.7" "1." "5.")
cut_pi_pt_lxy0to1_SS=("0.7" "1." "6.")
cut_pi_pt_lxy1to5_OS=("0.7" "2" "4.5")
cut_pi_pt_lxy1to5_SS=("0.7" "2" "4")
cut_pi_pt_lxygt5_OS=("2" "3." "4.")
cut_pi_pt_lxygt5_SS=("2" "3." "4.")

cut_sv_lxysig_lxy0to1_OS=("50" "90" "30")
cut_sv_lxysig_lxy0to1_SS=("50" "90" "30")
cut_sv_lxysig_lxy1to5_OS=("200" "270" "75")
cut_sv_lxysig_lxy1to5_SS=("150" "200" "75")
cut_sv_lxysig_lxygt5_OS=("230" "250" "150")
cut_sv_lxysig_lxygt5_SS=("230" "200" "150")

cut_min_mu_pi_dxysig_lxy0to1_OS=("0" "40" "0")
cut_min_mu_pi_dxysig_lxy0to1_SS=("0" "40" "0")
cut_min_mu_pi_dxysig_lxy1to5_OS=("20" "25" "0")
cut_min_mu_pi_dxysig_lxy1to5_SS=("20" "25" "0")
cut_min_mu_pi_dxysig_lxygt5_OS=("60" "40" "0")
cut_min_mu_pi_dxysig_lxygt5_SS=("60" "20" "0")

cut_hnl_cos2d_lxy0to1_OS=("1e-3" "1e-3" "1")
cut_hnl_cos2d_lxy0to1_SS=("1e-3" "1e-3" "1")
cut_hnl_cos2d_lxy1to5_OS=("1e-4" "1e-4" "1")
cut_hnl_cos2d_lxy1to5_SS=("1e-4" "1e-4" "1")
cut_hnl_cos2d_lxygt5_OS=("1e-5" "2e-5" "1")
cut_hnl_cos2d_lxygt5_SS=("1e-5" "2e-5" "1")

cut_mu_dxysig_lxy0to1_OS=("5" "15" "0")
cut_mu_dxysig_lxy0to1_SS=("5" "15" "0")
cut_mu_dxysig_lxy1to5_OS=("10" "20" "0")
cut_mu_dxysig_lxy1to5_SS=("10" "20" "0")
cut_mu_dxysig_lxygt5_OS=("25" "40" "0")
cut_mu_dxysig_lxygt5_SS=("25" "40" "0")

cut_pi_dxysig_lxy0to1_OS=("0" "20" "30")
cut_pi_dxysig_lxy0to1_SS=("0" "20" "30")
cut_pi_dxysig_lxy1to5_OS=("0" "0" "0")
cut_pi_dxysig_lxy1to5_SS=("0" "0" "0")
cut_pi_dxysig_lxygt5_OS=("45" "0" "0")
cut_pi_dxysig_lxygt5_SS=("45" "0" "0")

cut_max_mu_pi_dxysig_lxy0to1_OS=("0" "40" "20")
cut_max_mu_pi_dxysig_lxy0to1_SS=("0" "40" "20")
cut_max_mu_pi_dxysig_lxy1to5_OS=("40" "25" "10")
cut_max_mu_pi_dxysig_lxy1to5_SS=("40" "25" "10")
cut_max_mu_pi_dxysig_lxygt5_OS=("45" "20" "25")
cut_max_mu_pi_dxysig_lxygt5_SS=("45" "20" "25")


for MASS in "1" "3" "4p5"
do
  outdirlabel=$dirlabel'_m'${MASS}
  for CAT in "lxy0to1_OS" "lxy0to1_SS" "lxy1to5_OS" "lxy1to5_SS" "lxygt5_OS" "lxygt5_SS"
  do
    if [ ${MASS} == "1" ] ; then
      idx=0
    elif [ ${MASS} == "3" ] ; then
      idx=1
    elif [ ${MASS} == "4p5" ] ; then
      idx=2
    fi

    if [ ${CAT} == "lxy0to1_OS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxy0to1_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy0to1_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy0to1_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy0to1_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy0to1_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy0to1_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy0to1_OS[$idx]}
    elif [ ${CAT} == "lxy0to1_SS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxy0to1_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy0to1_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy0to1_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy0to1_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy0to1_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy0to1_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy0to1_SS[$idx]}
    elif [ ${CAT} == "lxy1to5_OS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxy1to5_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy1to5_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy1to5_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy1to5_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy1to5_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy1to5_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy1to5_OS[$idx]}
    elif [ ${CAT} == "lxy1to5_SS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxy1to5_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy1to5_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy1to5_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy1to5_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy1to5_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy1to5_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy1to5_SS[$idx]}
    elif [ ${CAT} == "lxygt5_OS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxygt5_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxygt5_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxygt5_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxygt5_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxygt5_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxygt5_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxygt5_OS[$idx]}
    elif [ ${CAT} == "lxygt5_SS" ] ; then
      cut_pi_pt=${cut_pi_pt_lxygt5_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxygt5_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxygt5_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxygt5_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxygt5_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxygt5_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxygt5_SS[$idx]}
    fi

    #sbatch -p standard --account t3 -o ./logs/log_${MASS}_${CAT}.txt -e ./logs/log_${MASS}_${CAT}.txt --job-name=selection_${MASS}_${CAT} run_selection.sh $outdirlabel ${MASS} ${CAT} $quantity $action $cut_pi_pt $cut_sv_lxysig $cut_mu_dxysig $cut_pi_dxysig $cut_max_mu_pi_dxysig $cut_min_mu_pi_dxysig $cut_hnl_cos2d 
    sh run_selection.sh $outdirlabel ${MASS} ${CAT} $quantity $action $cut_pi_pt $cut_sv_lxysig $cut_mu_dxysig $cut_pi_dxysig $cut_max_mu_pi_dxysig $cut_min_mu_pi_dxysig $cut_hnl_cos2d 

  done
done



