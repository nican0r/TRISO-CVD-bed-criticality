# TRISO-CVD-bed-criticality job container for AWS Batch.
#
# Base: openmc/openmc:latest — official OpenMC image with the CLI and Python
# bindings preinstalled against HDF5. Nuclear data is NOT baked in; it is
# mounted from EFS at /xs (see aws/cloudformation.yaml).
FROM openmc/openmc:latest

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        awscli \
        pyyaml \
        numpy \
        pandas \
        boto3 \
        uncertainties \
        h5py

WORKDIR /work

COPY furnace/ /work/furnace/
COPY scripts/ /work/scripts/
COPY params.yaml /work/params.yaml
COPY cases/.triso_cache/ /work/cases/.triso_cache/

ENV PYTHONPATH=/work \
    OPENMC_CROSS_SECTIONS=/xs/endfb80_hdf5/cross_sections.xml

ENTRYPOINT ["python", "/work/scripts/batch/entrypoint.py"]
