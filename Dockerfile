# adapted from instructions at https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.html#preparing-your-dockerfile

FROM gitlab-registry.cern.ch/atlas/statanalysis:0-3

# the base image provides a user called "atlas" that has the correct uid (1000) for binder

RUN source /release_setup.sh && python3 -m pip install --no-cache-dir notebook jupyterlab
# will have installed things in /home/atlas/.local/bin ... need to ensure this is in the path
ENV PATH="/home/atlas/.local/bin:${PATH}"
# Make sure the contents of our repo are in ${HOME}
WORKDIR /home/atlas
COPY --chown=atlas . ${HOME}