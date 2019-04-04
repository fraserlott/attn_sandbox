# attn_sandbox
Fraser C Lott, Met Office Hadley Centre.
Crown Copyright 2014-19.
Event attribution sandbox created for the EUPHEME project,
incorporating code from the previous EUCLEIA project
first used at the Extreme Events Summer School, ICTP Trieste 2014.

To run:
bokeh serve eupheme_sandbox.py

Variables you will need to localise in that file:
linuxpath (set to a web-readable directory with write permissions on your file system);
webpath (the http path corresponding to the directory above);
modelroot (location of the NetCDFs, '/project/detn/' on the CDN or '/s3/informatics-eupheme/' on AWS.

This is currently the version for running locally. 
The Kubernetes call from the Jupyter notebook version will need to be added to this for full AWS functionality.
