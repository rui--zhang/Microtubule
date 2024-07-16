#!/bin/csh -f
#
#   Control script to submit multiple jobs on a cluster using
#   the Sun Grid Engine. Each job processes N particles. N is
#   specified in the 'mparameters' file as 'increment'.

cp mparameters mparameters_run

set start 	= `grep start_process mparameters_run | awk '{print $2}'`
set end   	= `grep end_process mparameters_run | awk '{print $2}'`
set first 	= `grep first_particle mparameters_run | awk '{print $2}'`
set last  	= `grep last_particle mparameters_run | awk '{print $2}'`
set incr  	= `grep increment mparameters_run | awk '{print $2}'`
set data_input 	= `grep data_input mparameters_run | awk '{print $2}'`
set apix	= `grep pixel_size mparameters_run | awk '{print $2}'`
set twist       = `grep TWIST mparameters_run | awk '{print $2}'`
set rise        = `grep RISE mparameters_run | awk '{print $2}'`
set ncls	= 1

set working_directory = `pwd`
set SCRATCH = `pwd`

mainloop:

mkdir junk
cp mparameters mparameters_run_${start}
mv ${data_input}_${start}_r${ncls}.par junk
mv ${data_input}_${start}_r${ncls}.shft junk

@ prev = $start - 1

set firstn = $first
@ lastn = $first + $incr - 1

goto alignment_done

alignment_done:

  set firstn = $first
  @ lastn = $first + $incr - 1

checkdone:

    sleep 5
    while ( $firstn <= $last )
	
      grep --binary-files=text "overall score" $SCRATCH/${data_input}_${start}_r${ncls}.par_${firstn}_$lastn >& /dev/null
      if ($status) goto checkdone

      echo "Particles $firstn to $lastn, finished....  "`date`
      if ($firstn == $first ) head -60 $SCRATCH/${data_input}_${start}_r${ncls}.par_${firstn}_${lastn} | grep --binary-files=text C >> ${working_directory}/${data_input}_${start}_r${ncls}.par

      grep -v C --binary-files=text $SCRATCH/${data_input}_${start}_r${ncls}.par_${firstn}_${lastn} >> ${working_directory}/${data_input}_${start}_r${ncls}.par
      grep -v C --binary-files=text $SCRATCH/${data_input}_${start}_r${ncls}.shft_${firstn}_${lastn} >> ${working_directory}/${data_input}_${start}_r${ncls}.shft
      #\rm $SCRATCH/${data_input}_${start}.par_${firstn}_${lastn} >& /dev/null

      @ firstn = $firstn + $incr
      @ lastn = $lastn + $incr
      if ( $lastn >= $last ) set lastn = $last
    end

cp ${data_input}_mult_refine_n_r${ncls}.log_1_* ${data_input}_${start}_mult_refine_save 
reconstruct:
echo "Calculating 3D structure...."
  cp ${data_input}_${start}_r${ncls}.par ${data_input}_${start}_r${ncls}.par-bak
  cp ${data_input}_${start}_r${ncls}.shft ${data_input}_${start}_r${ncls}.shft-bak
  \rm $SCRATCH/${data_input}_mult_reconstruct.log >& /dev/null
  \rm $SCRATCH/${data_input}_${start}_r${ncls}.shft_* >& /dev/null
  \rm $SCRATCH/${data_input}_${start}_r${ncls}.par_* >&/dev/null
  #\rm $SCRATCH/${data_input}_${start}_r${ncls}.mrc_* >&/dev/null
  \rm $SCRATCH/${data_input}_mult_refine_n_r${ncls}.log_* >& /dev/null
