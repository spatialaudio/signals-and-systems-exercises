# signals-and-systems-exercises

## Tutorial

**Continuous- and Discrete-Time Signals and Systems - A Tutorial with Computational Examples**

- https://github.com/spatialaudio/signals-and-systems-exercises

Lecture and tutorial are designed for International Standard Classification of Education (ISCED) level 6.

## Lecture

This Jupyter notebook based tutorial using Python is accompanying the lecture

**Continuous-Time and Discrete-Time Signals and Systems - Theory and Computational Examples**

- https://github.com/spatialaudio/signals-and-systems-lecture

Lecture and tutorial are designed for International Standard Classification of Education (ISCED) level 6.


## German Version

We are in the design process of a very detailed 12 units tutorial to support remote learning for our students within COVID-19 era, thus the initial version is in German. Translations to English are scheduled ASAP.

Please see the LaTex main file `tutorial_latex_deu/sig_sys_ex.tex`.
There are several graphics included, which are created by the provided Jupyter notebooks. These have an ID with hex numbers in file name so that notebooks can be linked to the exercises within the tex project.

We might wish to compile all notebooks at once, then we can use:

`python3 -m nbconvert --ExecutePreprocessor.kernel_name="mydsp" --execute --inplace *.ipynb **/*.ipynb **/**/*.ipynb`

## License

- Creative Commons Attribution 4.0 International License (CC BY 4.0) for text/graphics
- MIT License for software

## Referencing

Please cite this open educational resource (OER) project as

Frank Schultz,
*Continuous- and Discrete-Time Signals and Systems - A Tutorial with Computational Examples*,
University of Rostock,
https://github.com/spatialaudio/signals-and-systems-exercises,
Either: short SHA of the master branch commit to be cited. Year of commit.
Or: zenodo.org DOI. release year of DOI. (TBD)

## Versions / Tags / Branches

- summer term 2020 version: https://github.com/spatialaudio/signals-and-systems-exercises/releases/tag/v0.1
- we use the branch `master` for developing
- we use the branch `output` to provide notebooks with rendered output. **Note**, that we will hard-reset this branch from time to time when updates are being made. So, please don't rely on commits from this branch!

## Authorship

University of Rostock:

- Frank Schultz (concept, design, main author)
- Till Rettberg (concept, design)
- Sascha Spors (proof reading, concept)
- Matthias Geier (proof reading, technical advisor)
- Vera Erbes (proof reading)
