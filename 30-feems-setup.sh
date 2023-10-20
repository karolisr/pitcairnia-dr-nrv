#!/usr/bin/env bash

deactivate &> /dev/null
conda deactivate &> /dev/null

# macOS
# brew install gdal geos gsl suite-sparse proj
# export CFLAGS="-I/opt/homebrew/include"

# Ubuntu
# sudo apt install python3.10-venv gdal-bin gdal-data libgdal-dev python3-gdal libgdal-grass geos-bin libgeos-dev gsl-bin libgsl-dev libsuitesparse-dev proj-bin libproj-dev

TEMPDIR="${HOME}/feems-temp-dir"
VENVDIR="${HOME}/.venvs"

mkdir -p "${VENVDIR}"
mkdir -p "${TEMPDIR}"
cd "${TEMPDIR}"
rm -rf ./feems-repo

python3.10 -m venv "${VENVDIR}/venv-py310-feems"
source "${VENVDIR}/venv-py310-feems/bin/activate"

pip install --upgrade pip
pip install wheel
pip install ipython
pip install shapely --force-reinstall --no-cache-dir --no-binary shapely

# ----------------------------------------------------------------------------
tee requirements.txt <<EOF
setuptools
flake8
pep8
pytest
fiona
click
numpy
scipy
pyproj
scikit-learn
pandas-plink
msprime
statsmodels
openpyxl
xlrd
PyYAML
scikit-sparse
cartopy
matplotlib==3.4.2
networkx==2.4.0
shapely
EOF
# ----------------------------------------------------------------------------

pip install -r requirements.txt

git clone https://github.com/NovembreLab/feems feems-repo
mv feems-repo/requirements.txt feems-repo/requirements.txt.bak
cp requirements.txt feems-repo/requirements.txt
pip install ./feems-repo

pip install numpy==1.23.5

deactivate

cd

rm -rf "${TEMPDIR}"
