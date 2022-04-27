#!/bin/bash

# ${1} dirlabel
# ${2} mass
# ${3} category
# ${4} quantity
# ${5} action
# ${6} cut_pi_pt
# ${7} cut_sv_lxysig
# ${8} cut_mu_dxysig
# ${9} cut_pi_dxysig
# ${10} cut_max_mu_pi_dxysig
# ${11} cut_min_mu_pi_dxysig
# ${12} cut_hnl_cos2d

homedir=$PWD

workdir='/scratch/anlyon/selection/'${1}'/'${2}'_'${3}'_'${5}
echo $workdir

mkdir -p $workdir
cp selection.py $workdir
cp tools.py $workdir
cp compute_yields.py $workdir
cp decays.py $workdir
cp decay_tools.py $workdir
cp -r ../objects $workdir/..
cp -r ../data $workdir/..

cd $workdir

echo "running script"
python selection.py --dirlabel ${1} --mass ${2} --category ${3} --quantity ${4} --action ${5} --cut_pi_pt ${6} --cut_sv_lxysig ${7} --cut_mu_dxysig ${8} --cut_pi_dxysig ${9} --cut_max_mu_pi_dxysig ${10} --cut_min_mu_pi_dxysig ${11} --cut_hnl_cos2d ${12}

echo "coyping the files"
cp -r ./myPlots/selection/${1} $homedir/myPlots/selection
echo "files copied"

cd $homedir
#rm -r $workdir/../objects
#rm -r $workdir/../data
rm -r $workdir

echo "End"
