# signals-and-systems-exercises

**Continuous- and Discrete-Time Signals and Systems - A Tutorial with Computational Examples**

This tutorial accompanies the lecture [Continuous-Time and Discrete-Time Signals and Systems - Theory and Computational Examples](https://github.com/spatialaudio/signals-and-systems-lecture). The lecture and the tutorial are designed for International Standard Classification of Education (ISCED) level 6 being held at University of Rostock, Germany.

Juypter notebooks can be accessed with the webserver services

- **dynamic** version using **mybinder**: https://mybinder.org/v2/gh/spatialaudio/signals-and-systems-exercises/output?filepath=index.ipynb
- **static** version using **nbviewer**: https://nbviewer.jupyter.org/github/spatialaudio/signals-and-systems-exercises/blob/output/index.ipynb

All sources (tex, ipynb, ...) can be found at github
- https://github.com/spatialaudio/signals-and-systems-exercises

## Anaconda Environment

The [Anaconda distribution](https://www.anaconda.com/distribution/) is a very convenient solution to install a required environment, i.e. to have access to the Jupyter Notebook renderer with a Python interpreter on a personal computer. It is very likely that a very recent installation of Anaconda already delivers all required packages just using the `base` environment. We actually do not need very special packages. However, if `base` is not working immediately, creating and activating a dedicated environment `mydsp` might be useful. Do the following steps

- get at least python 3.9x, numpy, sympy, scipy, matplotlib, notebook, jupyterlab, ipykernel, the other packages are very useful tools for convenience

`conda create -n mydsp python=3.9 pip numpy sympy scipy matplotlib notebook ipykernel jupyterlab pydocstyle pycodestyle autopep8 flake8 nb_conda jupyter_nbextensions_configurator jupyter_contrib_nbextensions`

- activate this environment with

`conda activate mydsp`

- Jupyter notebook renderer needs to know that we want to use the new Python environment

`python -m ipykernel install --user --name mydsp --display-name "mydsp"`

- get into the folder where the exercises are located, e.g.

`cd my_signals_and_systems_exercises_folder`

- start either a Jupyter notebook or Jupyter lab working environment via a local server instance by either

`jupyter notebook` or `jupyter lab`

If the above steps still lead to problems, the following lines created the working environment `mydsp`
- using `conda 4.12.0`, `conda-build 3.21.8`
- `conda create -n mydsp python=3.9.12 pip=22.0.4 numpy=1.22.3 sympy=1.10.1 scipy=1.8.0 matplotlib=3.5.1 notebook=6.4.10 ipykernel=6.11.0 jupyterlab=3.3.2 pydocstyle=6.1.1 pycodestyle=2.8.0 autopep8=1.6.0 flake8=4.0.1 nb_conda=2.2.1 jupyter_nbextensions_configurator=0.4.1 jupyter_contrib_nbextensions=0.5.1`
- `pip install soundfile`

## German Version

We are still in the design and polishing process of a very detailed 12 units tutorial to support remote learning for our students within COVID-19 era, thus the initial version is in German. Translations to English are scheduled ASAP. Please see the LaTex main file `tutorial_latex_deu/sig_sys_ex.tex`. There are several graphics included, which are created by the provided Jupyter notebooks. These have an ID with hex numbers in file name, so that notebooks can be linked to the exercises within the tex project. We might wish to compile all notebooks at once, then we can use
`python3 -m nbconvert --ExecutePreprocessor.kernel_name="mydsp" --execute --inplace *.ipynb **/*.ipynb **/**/*.ipynb`

## License

- Creative Commons Attribution 4.0 International License (CC BY 4.0) for text/graphics
- MIT License for software

## Versions / Tags / Branches

- summer term 2021 version: https://github.com/spatialaudio/signals-and-systems-exercises/releases/tag/v0.2
- summer term 2020 version: https://github.com/spatialaudio/signals-and-systems-exercises/releases/tag/v0.1
- we use the branch `master` for developing
- we use the branch `output` to provide notebooks with rendered output. **Note**, that we will hard-reset this branch from time to time when updates are being made. So, please do not rely on commits from this branch!

## Referencing

Please cite this open educational resource (OER) project as
*Frank Schultz, Continuous- and Discrete-Time Signals and Systems - A Tutorial Featuring Computational Examples, University of Rostock* with ``github URL, commit number and/or version tag, year, (dedicated files and/or content)``.

## Authorship

University of Rostock:

- Frank Schultz (concept, design, main author)
- Robert Hauser (proof reading, concept, author)
- Sascha Spors (proof reading, concept)
- Till Rettberg (concept, design)
- Matthias Geier (proof reading, technical advisor)
- Vera Erbes (proof reading)
