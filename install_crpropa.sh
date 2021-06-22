# Install virtualenv to create a python virtual environment
sudo pip3 install virtualenv -y

# Create a virtual environment with python 2.7 and activate it
export CRPROPA_DIR=$HOME"/.virtualenvs/crpropa"
mkdir -p $CRPROPA_DIR
virtualenv --python=/usr/bin/python2.7 $CRPROPA_DIR
source $CRPROPA_DIR"/bin/activate"
cd $CRPROPA_DIR

# Install required dependencies
pip install matplotlib numpy pandas
sudo apt install build-essential git cmake swig gfortran python-dev fftw3-dev zlib1g-dev libmuparser-dev libhdf5-dev pkg-config

# Download crpropa from repository and install it
git clone https://github.com/CRPropa/CRPropa3.git
cd CRPropa3
mkdir build
cd build
sudo CMAKE_PREFIX_PATH=$CRPROPA_DIR cmake -DCMAKE_INSTALL_PREFIX=$CRPROPA_DIR ..
sudo make
sudo make test
sudo make install
