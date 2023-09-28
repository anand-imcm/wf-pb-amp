FROM mambaorg/micromamba:1.4.2

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/*_requirements.* /tmp/

COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts /scripts

RUN micromamba install -y -n base -f /tmp/conda_requirements.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN pip install -r /tmp/python_requirements.txt

RUN python /scripts/fix_py_vcf_parser.py

# adding the path to the tools managed by mamba
ENV PATH="/opt/conda/bin:$PATH"