# ASTEC

This file contains instructions for installing and running the ASTEC reconstruction algorithm.

ASTEC is a **Segmentation and Tracking algorithm from Contact-dependent cell communications drive morphological invariance during ascidian embryogenesis**.

This work was originaly published on **bioRxiv** in 2018:

**Contact-dependent cell-cell communications drive morphological invariance during ascidian embryogenesis**, _Léo Guignard, Ulla-Maj Fiuza, Bruno Leggio, Emmanuel Faure, Julien Laussu, Lars Hufnagel, Grégoire Malandain, Christophe Godin, Patrick Lemaire_, bioRxiv 2018; doi: https://doi.org/10.1101/238741.

It was later published in **Science** in 2020:

**Contact area–dependent cell communication and the morphological invariance of ascidian embryogenesis**, _Léo Guignard, Ulla-Maj Fiuza, Bruno Leggio, Julien Laussu, Emmanuel Faure, Gaël Michelin, Kilian Biasuz, Lars Hufnagel, Grégoire Malandain, Christophe Godin, Patrick Lemaire_, Science 2020; doi: https://doi.org/10.1126/science.aar5663.


## Installation

Only works for Linux or MacOs systems.

### I - User installation (without git)

Requires `conda`.

1. Create a conda environment (here the environment is named
`astec-test`)

	```bash
	conda create -n astec -c morpheme -c conda-forge astec
	```

2. Activate the built conda environment

	```bash
	conda activate astec
	```

### II - User installation (with git)

Requires `conda` and `git`.

1. Astec code can be found at
   [gitlab.inria.fr/astec/astec](http://gitlab.inria.fr/astec/astec). It
   can be downloaded with

	```bash
	git clone https://gitlab.inria.fr/astec/astec.git
	```

	It creates an `astec` directory.

2. Create a conda environment named `astec`

	```bash
	cd astec
	conda env create -f pkg/env/astec.yaml
	```

3. Activate the built conda environment

	```bash
	conda activate astec
	```

### III - Developer  installation (with git)

Requires `conda` and `git`.

1. Astec code can be found at
   [gitlab.inria.fr/astec/astec](http://gitlab.inria.fr/astec/astec). It
   can be downloaded with

	```bash
	git clone https://gitlab.inria.fr/astec/astec.git
	```

	It creates an `astec` directory.

2. Create a conda environment named `astec`

	```bash
	cd astec
	conda env create -f pkg/env/astec-dev.yaml
	```

3. Activate the built conda environment

	```bash
	conda activate astec-dev
	```
	
4. Install astec package for use

	```bash
	python -m pip install -e .
	```

	The `-e` option install the package in "editable" mode, this is
    want you want if you aim at contributing to the astec
    project. This last command has to be repeated (within the conda
    environment every time the astec code has been modified).


## Tutorial

A tutorial can be found at [gitlab.inria.fr/astec/astec-tutorial](https://gitlab.inria.fr/astec/astec-tutorial).

## Documentation

A documentation can be found at [astec.gitlabpages.inria.fr/astec/](https://astec.gitlabpages.inria.fr/astec/).
