FROM climateimpactlab/miniconda3-gfortran-netcdf

WORKDIR /

RUN apt-get install --yes libproj-dev proj-data proj-bin;
RUN apt-get install --yes libgeos-dev;
RUN pip install cython;
RUN pip install \
        "numpy>=1.6" \
        "matplotlib>=1.0.1" \
        "scipy>=0.10.0" \
        "pandas>=0.24.2" \
        "xarray>=0.12.0" \
        "netcdf4>=1.5.0.1" \
        "nose==1.3.7" \
        "geopy>=1.19.0";

RUN git clone --branch=master --depth 100 --quiet git://github.com/ClimateImpactLab/clawpack
WORKDIR /clawpack
RUN git submodule init
RUN git submodule update -j 6

RUN pip install -e .

ENV PYTHONPATH /clawpack:$PYTHONPATH
ENV CLAW /clawpack
ENV FC gfortran

WORKDIR /home
RUN mkdir /home/examples

CMD "/bin/bash"
