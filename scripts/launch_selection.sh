#!/bin/bash

#dirlabel='selection_significance_diff_dataD1_v3'
dirlabel='study_central_V12_08Aug22_10files'
quantity='b_mass'
action='scan_significance_diff'
submit_batch='1'

# list of cuts
# syntax (m1 m3 m4p5)
cut_b_mass_lxy0to1_OS=("0.7" "0.7" "0.7")
cut_b_mass_lxy0to1_SS=("0.7" "0.7" "0.7")
cut_b_mass_lxy1to5_OS=("0.7" "0.7" "0.7")
cut_b_mass_lxy1to5_SS=("0.7" "0.7" "0.7")
cut_b_mass_lxygt5_OS=("0.7" "0.7" "0.7")
cut_b_mass_lxygt5_SS=("0.7" "0.7" "0.7")

cut_pi_pt_lxy0to1_OS=("0.7" "0.8" "2.")
cut_pi_pt_lxy0to1_SS=("0.7" "0.8" "2.")
cut_pi_pt_lxy1to5_OS=("0.7" "1.5" "3.")
cut_pi_pt_lxy1to5_SS=("0.7" "1.5" "3.")
cut_pi_pt_lxygt5_OS=("1.5" "2.5" "4.")
cut_pi_pt_lxygt5_SS=("1.5" "2.5" "4.")

cut_sv_lxysig_lxy0to1_OS=("30" "90" "100")
cut_sv_lxysig_lxy0to1_SS=("50" "100" "70")
cut_sv_lxysig_lxy1to5_OS=("150" "300" "100")
cut_sv_lxysig_lxy1to5_SS=("150" "300" "100")
cut_sv_lxysig_lxygt5_OS=("150" "300" "200")
cut_sv_lxysig_lxygt5_SS=("150" "300" "200")

cut_min_mu_pi_dxysig_lxy0to1_OS=("0" "5" "0")
cut_min_mu_pi_dxysig_lxy0to1_SS=("0" "5" "0")
cut_min_mu_pi_dxysig_lxy1to5_OS=("10" "60" "0")
cut_min_mu_pi_dxysig_lxy1to5_SS=("10" "60" "0")
cut_min_mu_pi_dxysig_lxygt5_OS=("30" "80" "0")
cut_min_mu_pi_dxysig_lxygt5_SS=("30" "80" "0")

cut_hnl_cos2d_lxy0to1_OS=("2e-3" "2e-3" "1")
cut_hnl_cos2d_lxy0to1_SS=("2e-3" "2e-3" "1")
cut_hnl_cos2d_lxy1to5_OS=("1e-4" "2e-4" "1")
cut_hnl_cos2d_lxy1to5_SS=("3e-4" "2e-4" "1")
cut_hnl_cos2d_lxygt5_OS=("3e-5" "2e-5" "1")
cut_hnl_cos2d_lxygt5_SS=("3e-5" "2e-5" "1")

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
      cut_b_mass=${cut_b_mass_lxy0to1_OS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxy0to1_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy0to1_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy0to1_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy0to1_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy0to1_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy0to1_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy0to1_OS[$idx]}
    elif [ ${CAT} == "lxy0to1_SS" ] ; then
      cut_b_mass=${cut_b_mass_lxy0to1_SS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxy0to1_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy0to1_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy0to1_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy0to1_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy0to1_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy0to1_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy0to1_SS[$idx]}
    elif [ ${CAT} == "lxy1to5_OS" ] ; then
      cut_b_mass=${cut_b_mass_lxy1to5_OS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxy1to5_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy1to5_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy1to5_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy1to5_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy1to5_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy1to5_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy1to5_OS[$idx]}
    elif [ ${CAT} == "lxy1to5_SS" ] ; then
      cut_b_mass=${cut_b_mass_lxy1to5_SS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxy1to5_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxy1to5_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxy1to5_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxy1to5_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxy1to5_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxy1to5_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxy1to5_SS[$idx]}
    elif [ ${CAT} == "lxygt5_OS" ] ; then
      cut_b_mass=${cut_b_mass_lxygt5_OS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxygt5_OS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxygt5_OS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxygt5_OS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxygt5_OS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxygt5_OS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxygt5_OS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxygt5_OS[$idx]}
    elif [ ${CAT} == "lxygt5_SS" ] ; then
      cut_b_mass=${cut_b_mass_lxygt5_SS[$idx]}
      cut_pi_pt=${cut_pi_pt_lxygt5_SS[$idx]}
      cut_sv_lxysig=${cut_sv_lxysig_lxygt5_SS[$idx]}
      cut_mu_dxysig=${cut_mu_dxysig_lxygt5_SS[$idx]}
      cut_pi_dxysig=${cut_pi_dxysig_lxygt5_SS[$idx]}
      cut_max_mu_pi_dxysig=${cut_max_mu_pi_dxysig_lxygt5_SS[$idx]}
      cut_min_mu_pi_dxysig=${cut_min_mu_pi_dxysig_lxygt5_SS[$idx]}
      cut_hnl_cos2d=${cut_hnl_cos2d_lxygt5_SS[$idx]}
    fi

    if [ $submit_batch == '1' ] ; then
      sbatch -p standard --account t3 -o ./logs/log_${MASS}_${CAT}.txt -e ./logs/log_${MASS}_${CAT}.txt --job-name=selection_${MASS}_${CAT} --time 02:00:00 run_selection.sh $outdirlabel ${MASS} ${CAT} $quantity $action $cut_b_mass $cut_pi_pt $cut_sv_lxysig $cut_mu_dxysig $cut_pi_dxysig $cut_max_mu_pi_dxysig $cut_min_mu_pi_dxysig $cut_hnl_cos2d 
    elif [ $submit_batch == '0' ] ; then
      sh run_selection.sh $outdirlabel ${MASS} ${CAT} $quantity $action $cut_b_mass $cut_pi_pt $cut_sv_lxysig $cut_mu_dxysig $cut_pi_dxysig $cut_max_mu_pi_dxysig $cut_min_mu_pi_dxysig $cut_hnl_cos2d 
    fi

  done
done



