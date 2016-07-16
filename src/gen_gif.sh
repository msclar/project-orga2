#!/bin/bash
	# receive parameters
	max_i=$1
	max_j=$2
	delta_t=$3
	max_cant_delta_t=$4
	cant_prints=$5
	gif_name=$6

	# remove folders with old data
	rm -r output
	mkdir -p output
	
	rm -r images
	mkdir -p images
	
	# execute jacobi, plot with matlab and create gif
	./jacobi-optimizado $max_i $max_j $delta_t $max_cant_delta_t $cant_prints
	#alias matlab='/usr/local/MATLAB/R2016a/bin/matlab'
	/usr/local/MATLAB/R2016a/bin/matlab -nodisplay -r "try, graficar, catch, end, quit"
	convert -delay 20 -loop 0 $(ls -1v images/*.png) $gif_name.gif
