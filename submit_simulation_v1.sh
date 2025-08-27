python3 geomConvert.py 

docker run --rm --tty --interactive --volume "${PWD}:/wd" --env OMP_NUM_THREADS={n_threads} eisenforschung/damask-grid:latest --load tensionX.yaml --geom grainsID.vti
python3 postProcessXX.py

docker run --rm --tty --interactive --volume "${PWD}:/wd" --env OMP_NUM_THREADS={n_threads} eisenforschung/damask-grid:latest --load tensionY.yaml --geom grainsID.vti
python3 postProcessYY.py


docker run --rm --tty --interactive --volume "${PWD}:/wd" --env OMP_NUM_THREADS={n_threads} eisenforschung/damask-grid:latest --load tensionZ.yaml --geom grainsID.vti
python3 postProcessZZ.py

