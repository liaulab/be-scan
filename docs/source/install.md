# Installation

## Installation

To install BE-SCAN:

```console
$ pip install git+https://github.com/liaulab/be-scan.git@simple_version
```

## Conda Environment

To setup an environment with BE-SCAN: 

Create a conda environment with: 
```console
$ conda create --name bescan-test python=3.9
```

Activate the environment with: 
```console
$ conda activate bescan-test 
```

Install be-scan with: 
```console
$ python -m pip install git+https://github.com/liaulab/be-scan.git@simple_version
``````

## Using Jupyter Notebook

To use BE-SCAN on jupyter notebook:

If you intent to use the jupyter notebook interface, install jupyter with: 
```console
$ conda install jupyterlab
```
```console
$ pip install jupyterlab
```

Add conda environment to jupyter notebook:
```console
$ python -m ipykernel install --user --name=bescan-test
```
Activate jupyter notebook with: 
```console
$ jupyter notebook
```
Select be-scan-test as kernel
