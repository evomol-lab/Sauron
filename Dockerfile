# Sauron: residue interaction network (RIN) calculator, web interface + command line.
#
# Image names are fully qualified (docker.io/..., ghcr.io/...) so the same file
# builds with Docker and with Podman, which does not resolve short names by default.
#
#   podman build -t docker.io/evomol/sauron:<version> .
#   podman run --rm -p 5000:5000 docker.io/evomol/sauron:<version>          # web interface
#   podman run --rm -v "$PWD":/data -w /data docker.io/evomol/sauron:<version> \
#       python /opt/sauron/sauron.py structure.cif                          # command line

ARG PYTHON_IMAGE=docker.io/library/python:3.12-slim

# ---------------------------------------------------------------- build stage
# requirements.txt pins the tested environment (uv.lock) and installs CPU-only
# PyTorch; the same file is used by Apptainer.def and for local installs.
FROM ${PYTHON_IMAGE} AS build

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN python -m venv /opt/venv

COPY requirements.txt /tmp/requirements.txt
RUN /opt/venv/bin/pip install -r /tmp/requirements.txt

# ---------------------------------------------------------------- runtime stage
FROM ${PYTHON_IMAGE}

LABEL org.opencontainers.image.title="Sauron" \
      org.opencontainers.image.description="Residue interaction network (RIN) calculator: web interface and command line" \
      org.opencontainers.image.source="https://github.com/jpmslima/Sauron" \
      org.opencontainers.image.vendor="EvoMol-Lab, BioMe, UFRN" \
      org.opencontainers.image.licenses="PolyForm-Noncommercial-1.0.0"

ENV PATH=/opt/venv/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    SAURON_UPLOAD_DIR=/opt/sauron/uploads

COPY --from=build /opt/venv /opt/venv

WORKDIR /opt/sauron
COPY sauron.py SauronGUI.py DOCUMENTATION.md README.md LICENSE ./
COPY templates/ templates/
COPY static/ static/

# Results of the web interface; mount a host directory here to keep them.
RUN mkdir -p uploads && chmod 777 uploads

EXPOSE 5000

# Python running as PID 1 ignores SIGTERM; SIGINT stops the Flask server cleanly.
STOPSIGNAL SIGINT

# Absolute path: Apptainer images built from this one (docker://...) ignore WORKDIR.
CMD ["python", "/opt/sauron/SauronGUI.py"]
