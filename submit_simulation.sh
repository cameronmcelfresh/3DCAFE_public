python geomConvert.py 

docker run --rm --tty --interactive --volume "${PWD}:/wd" --env OMP_NUM_THREADS={n_threads} eisenforschung/damask-grid:latest --load tensionX.yaml --geom grainsID.vti

python3 DAMASK_postProcess.py

rm grainsID_tensionX.hdf5
rm grainsID_tensionX.C_ref
rm grainsID_tensionX.sta