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


ENV PYTHONPATH /clawpack:$PYTHONPATH
ENV CLAW /clawpack
ENV FC gfortran
ENV NETCDF4_DIR /usr/local

WORKDIR /home
RUN mkdir /home/examples

COPY prepare.sh /usr/bin/prepare.sh
RUN chmod +x /usr/bin/prepare.sh

ENTRYPOINT ["tini", "--", "/usr/bin/prepare.sh"]
CMD "/bin/bash"
