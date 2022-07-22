#!/bin/bash

compile_files(){
	cd LinuxC++Compile
	g++ -fPIC -fpermissive -w -c NRpyDNAcode.cpp -o NRpyDNAcode.o -I/usr/include/python3.8 \
	 -I/usr/local/lib/python3.8/dist-packages/numpy/core/include
	g++ -shared NRpyDNAcode.o -o NRpyDNAcode.so

	g++ -fPIC -fpermissive -w -c NRpyRS.cpp -o NRpyRS.o -I/usr/include/python3.8 \
	 -I/usr/local/lib/python3.8/dist-packages/numpy/core/include
	g++ -shared NRpyRS.o -o NRpyRS.so
	mv NRpyDNAcode.so ../
	mv NRpyRS.so ../
	echo "compilation done"
	cd ..
}

run_internal(){
	# run python script with Environment Variable - PYTHONMALLOC=malloc
	PYTHONMALLOC=malloc python3.8 project.py
}

run(){
	# compile C++ files
	compile_files
	# run python script
	run_internal
}

