# hdag-benchmark

## Dependencies

Install dependencies via:

    conda env create -f environment.yml
    conda activate hdag-benchmark


## Install

    pip install -e .


## Additional Requirements

We also make use of beast2_xml, a python package for loading fasta files into XML format.
Unfortunately, this package is incompatible with protobuf.
So, we create a new conda environment to switch to when loading the XML files.

    conda create -n beast-xml
    conda activate beast-xml
    pip install beast2_xml