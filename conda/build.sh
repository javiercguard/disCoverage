echo "Installing disCoverage"
cp preparePlottingCoverage.py "$CONDA_PREFIX"/lib/python3.8/site-packages
mkdir -p "$CONDA_PREFIX"/supp
cp hg38.bed "$CONDA_PREFIX"/supp
cp disCoverage.R "$CONDA_PREFIX"/supp
