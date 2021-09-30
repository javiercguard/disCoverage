echo "Installing disCoverage"
$PYTHON setup.py install
cp hg38.bed $CONDA_PREFIX/
cp disCoverage.R $CONDA_PREFIX/
