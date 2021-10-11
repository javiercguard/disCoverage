echo "Installing disCoverage"
$PYTHON setup.py install
mkdir -p "$CONDA_PREFIX"/supp
cp hg38.bed "$CONDA_PREFIX"/supp
cp disCoverage.R "$CONDA_PREFIX"/supp
