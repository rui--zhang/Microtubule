#!/bin/csh -f
#unlimit
limit coredumpsize 0
set working_directory	= `pwd`
cp mparameters mparameters_run
set SCRATCH		= `pwd`
cd $SCRATCH
echo 'Did I add SSNR table to the parameter file?'

set bin_dir		= /home/ruiz/bin/frealign_v9.11/bin

set start = $3
set data_input	= `grep data_input mparameters_run | awk '{print $2}'`
set raw_images	= `grep raw_images_high mparameters_run | awk '{print $2}'`
if ( $raw_images == "" ) set raw_images = `grep raw_images mparameters_run | awk '{print $2}'`
set raw_images	= `echo ${raw_images:r}`
set extension	= `ls $raw_images.* | head -1`
if ( $extension == "" ) then
  set raw_images = ${working_directory}/${raw_images}
  set extension	= `ls $raw_images.* | head -1`
endif
set extension	= `echo ${extension:e}`
set mem_per_cpu	= `grep mem_per_cpu mparameters_run | awk '{print $2}'`
if ( $status || $mem_per_cpu == "" ) then
  set mem_per_cpu	= 4096
endif
set target	= `grep thresh_refine mparameters_run | awk '{print $2}'`
set pbc		= `grep PBC mparameters_run | awk '{print $2}'`
# set boff	= `grep BOFF mparameters_run | awk '{print $2}'`
set dang	= `grep DANG mparameters_run | awk '{print $2}'`
set itmax	= `grep ITMAX mparameters_run | awk '{print $2}'`
set mode	= `grep MODE mparameters_run | awk '{print $2}'`
set FMAG	= `grep FMAG mparameters_run | awk '{print $2}'`
set FDEF	= `grep FDEF mparameters_run | awk '{print $2}'`
set FASTIG	= `grep FASTIG mparameters_run | awk '{print $2}'`
set FPART	= `grep FPART mparameters_run | awk '{print $2}'`
set FMATCH	= `grep FMATCH mparameters_run | awk '{print $2}'`
set FBOOST	= `grep FBOOST mparameters_run | awk '{print $2}'`
set rref	= `grep res_high_refinement mparameters_run | awk '{print $2}'`
set rclas	= `grep res_high_class mparameters_run | awk '{print $2}'`
set rlowref	= `grep res_low_refinement mparameters_run | awk '{print $2}'`
set rbf		= `grep RBfactor mparameters_run | awk '{print $2}'`
set sym		= `grep Symmetry mparameters_run | awk '{print $2}'`
set alpha	= `grep ALPHA mparameters_run | awk '{print $2}'`
set rise	= `grep RISE mparameters_run | awk '{print $2}'`
set nsub	= `grep NSUBUNITS mparameters_run | awk '{print $2}'`
set nstart	= `grep NSTARTS mparameters_run | awk '{print $2}'`
set stiff	= `grep STIFFNESS mparameters_run | awk '{print $2}'`
set pix		= `grep pix_size mparameters_run | awk '{print $2}'`
set kV		= `grep Voltage mparameters_run | awk '{print $2}'`
set AmpC	= `grep Amp_contrast mparameters_run | awk '{print $2}'`
set ImC		= `grep image_contrast mparameters_run | awk '{print $2}'`
if ( ! $status ) then
  if ( $ImC == "P" ) set AmpC = `echo -$AmpC`
endif
set XSTD	= `grep XSTD mparameters_run | awk '{print $2}'`
set dstep	= `grep dstep mparameters_run | awk '{print $2}'`
set ro		= `grep outer_radius mparameters_run | awk '{print $2}'`
set ri		= `grep inner_radius mparameters_run | awk '{print $2}'`
set rlowref = `echo $rlowref $ro | awk '{if ($1+0.0 == 0.0) {print 2.0*$2} else {print $1} }'`
set MW		= `grep mol_mass mparameters_run | awk '{print $2}'`
set cs		= `grep Aberration mparameters_run | awk '{print $2}'`
#set mode	= ${5}
set psi		= ${6}
set theta	= ${7}
set phi		= ${8}
set shx		= ${9}
set shy		= ${10}
#
set form = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text Opening | awk '{print $2}'`
set fm = "M"
if ( $form == "SPIDER" ) set fm = "S"
if ( $form == "IMAGIC" ) set fm = "I"
#
set imem = 3
set nx = `${bin_dir}/fheader.exe ${raw_images}.${extension} | grep --binary-files=text NX | awk '{print $4}'`
set mem_big = `echo $nx | awk '{print int(10 * $1^3 * 4 * 66 /1024^3 + 1)/10}'`
if ( `echo $mem_big | awk '{print int(1024 * $1)}'` > $mem_per_cpu ) then
  set imem = 0
endif

if (${mode} == "") set mode = 1

### check for shift & angle refinement
if (${psi} == "") set psi = 1
if (${theta} == "") set theta = 1
if (${phi} == "") set phi = 1
if (${shx} == "") set shx = 1
if (${shy} == "") set shy = 1

set ifsc = 0
if ( $FBOOST == "T" ) set ifsc = -1

@ prev = $start - 1
#
\rm ${data_input}_${start}_r${4}.par_${1}_${2} >& /dev/null
#
time ${bin_dir}/frealign_v9.exe << eot >& ${data_input}_mult_refine_n_r${4}.log_${1}_${2}
M,${mode},${FMAG},${FDEF},${FASTIG},${FPART},0,F,F,F,${FMATCH},${ifsc},F,0,1	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
${ro},${ri},${pix},${MW},${AmpC},${XSTD},${pbc},0.0,${dang},${itmax},20	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
${psi},${theta},${phi},${shx},${shy}					!MASK
${1},${2}								!IFIRST,ILAST
0									!ASYM symmetry card (I=icosahedral)
1.0,${dstep},${target},0.0,${cs},${kV},0.0,0.0				!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rref},   ${rlowref},   ${rref}, ${rclas}, 100.0, ${rbf}		!RREC,RMIN,RMAX,RCLAS,DFSTD,RBFACT
${raw_images}.${extension}
${data_input}_reproject_r${4}.${extension}_${1}_${2}
${working_directory}/${data_input}_${prev}_r${4}.par
${data_input}_${start}_r${4}.par_${1}_${2}
${data_input}_${start}_r${4}.shft_${1}_${2}
-100., 0., 0., 0., 0., 0., 0., 0.					! terminator with RELMAG=0.0
${working_directory}/${data_input}_${prev}_r${4}.${extension}
${data_input}_${start}_r${4}_weights
${data_input}_${start}_r${4}_map1.${extension}
${data_input}_${start}_r${4}_map2.${extension}
${data_input}_${start}_r${4}_phasediffs
${data_input}_${start}_r${4}_pointspread
eot
#
#echo Job on $HOST finished >> stderr
#
