from distutils.core import setup

setup(
	name = "disCoverage",
	version = "1.0",
	py_modules = ['preparePlottingCoverage'],
	scripts = ['disCoverage'],
	package_data = {'disCoverage': ['hg38.bed', 'disCoverage.R']},
    including_package_data = True,
	)