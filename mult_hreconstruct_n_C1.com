#!/bin/csh -f
#unlimit
limit coredumpsize 0
set working_directory	= `pwd`
set SCRATCH		= `pwd`
cd $SCRATCH

set bin_dir		= /home/zhangrui/bin/frealign_v9.09/bin

set start = $3
cp mparameters_run mparameters_run_${start}
@ prev = $start - 1
set data_input	= `grep data_input mparameters_run | awk '{print $2}'`
set raw_images	= `grep raw_images_low mparameters_run | awk '{print $2}'`
if ( $raw_images == "" ) set raw_images = `grep raw_images mparameters_run | awk '{print $2}'`
set raw_images	= `echo ${raw_images:r}`
set extension	= `ls $raw_images.* | head -1`
if ( $extension == "" ) then
  set raw_images = ${working_directory}/${raw_images}
  set extension	= `ls $raw_images.* | head -1`
endif
set extension	= `echo ${extension:e}`
set thresh	= `grep thresh_reconst mparameters_run | awk '{print $2}'`
set pbc		= `grep PBC mparameters_run | awk '{print $2}'`
# set boff	= `grep BOFF mparameters_run | awk '{print $2}'`
set dang	= `grep DANG mparameters_run | awk '{print $2}'`
set itmax	= `grep ITMAX mparameters_run | awk '{print $2}'`
set mode	= `grep MODE mparameters_run | awk '{print $2}'`
set FFILT	= `grep FFILT mparameters_run | awk '{print $2}'`
set FBEAUT	= `grep FBEAUT mparameters_run | awk '{print $2}'`
set rrec	= `grep res_reconstruction mparameters_run | awk '{print $2}'`
set rref	= `grep res_high_refinement mparameters_run | awk '{print $2}'`
set rclas	= `grep res_high_class mparameters_run | awk '{print $2}'`
set rbf		= `grep RBfactor mparameters_run | awk '{print $2}'`
set sym		= `grep Symmetry mparameters_run | awk '{print $2}'`
set alpha	= `grep ALPHA mparameters_run | awk '{print $2}'`
set rise	= `grep RISE mparameters_run | awk '{print $2}'`
set nsub	= `grep NSUBUNITS mparameters_run | awk '{print $2}'`
set nstart	= `grep NSTARTS mparameters_run | awk '{print $2}'`
set stiff	= `grep STIFFNESS mparameters_run | awk '{print $2}'`
set pix		= `grep pix_size mparameters_run | awk '{print $2}'`
set rrec = `echo $rrec $pix | awk '{if ($1+0.0 == 0.0) {print 2.1*$2} else {print $1} }'`
set kV		= `grep Voltage mparameters_run | awk '{print $2}'`
set AmpC	= `grep Amp_contrast mparameters_run | awk '{print $2}'`
set ImC		= `grep image_contrast mparameters_run | awk '{print $2}'`
if ( ! $status ) then
  if ( $ImC == "P" ) set AmpC = `echo -$AmpC`
endif
set dstep	= `grep dstep mparameters_run | awk '{print $2}'`
set ro		= `grep outer_radius mparameters_run | awk '{print $2}'`
set ri		= `grep inner_radius mparameters_run | awk '{print $2}'`
set MW		= `grep mol_mass mparameters_run | awk '{print $2}'`
set cs		= `grep Aberration mparameters_run | awk '{print $2}'`
set mp_cpus	= `grep mp_cpus mparameters_run | awk '{print $2}'`
#
set m = ""
if ( ${mp_cpus} > 1 ) then
  set m = "_mp"
  setenv NCPUS ${mp_cpus}
  setenv OMP_NUM_THREADS ${mp_cpus}
endif

set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
set fm = "M"
if ( $form == "SPIDER" ) set fm = "S"
if ( $form == "IMAGIC" ) set fm = "I"
#
\rm ${data_input}_${start}_r${4}_n${1}.${extension}
\rm ${data_input}_${start}_r${4}_n${1}.res

time ${bin_dir}/frealign_v9_mp.exe << eot >& ${data_input}_mult_reconstruct_r${4}_n${1}.log
M,0,F,F,F,F,0,${FBEAUT},${FFILT},F,F,0,F,2,1		!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
${ro},${ri},${pix},${MW},${AmpC},0.0,${pbc},0.0,10.,1,10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
1 1 1 1 1						!MASK
${1},${2}						!IFIRST,ILAST 
0							!ASYM symmetry card (I=icosahedral)
1.,${dstep},60.0,${thresh},${cs},${kV},0.0,0.0		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rrec}, 200.0, ${rref}, ${rclas}, 100.0, ${rbf}	!RREC,RMIN,RMAX,RCLAS,DFSTD,RBFACT
${raw_images}.${extension}
/dev/null
${working_directory}/${data_input}_${start}_r${4}.par
${data_input}_${start}_r${4}.res
${data_input}_${start}_r${4}_n${1}.shft
0., 0., 0., 0., 0., 0., 0., 0.				! terminator with RELMAG=0.0
${data_input}_${start}_r${4}_n${1}.${extension}
${data_input}_${start}_r${4}_weights
${data_input}_${start}_r${4}_map1.${extension}
${data_input}_${start}_r${4}_map2.${extension}
${data_input}_${start}_r${4}_phasediffs
${data_input}_${start}_r${4}_pointspread
eot
#
\rm ${data_input}_${start}_r${4}_weights
\rm ${data_input}_${start}_r${4}_n${1}.shft
\rm ${data_input}_${start}_r${4}_phasediffs
\rm ${data_input}_${start}_r${4}_pointspread
#
echo 'mreconstruct.com finished' >> ${data_input}_mult_reconstruct_r${4}_n${1}.log
date
