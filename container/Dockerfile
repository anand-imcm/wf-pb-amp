FROM mambaorg/micromamba:1.5.1

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/environment.yml /tmp/

COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts /scripts

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN pip3 install "pandas<=2.1.1" "pysam<=0.22.0"

# adding the path to the tools managed by mamba
ENV PATH="/opt/conda/bin:$PATH"