# molspotter

## Disclaimer

This repository is based on the awesome [`silly_walks`](https://github.com/PatWalters/silly_walks) work of Patrick Walters.

## Summary

***Identify silly molecules.***

`molspotter` helps finding unusual/silly atom arrangements in molecules.
To do so a reference dataset of molecules is used. Any fingerprint bit of a new molecule not appearing in the reference set is considered silly. 

## Installation

```commandline
pip install molspotter
```

## Usage

Three prepared spotters are included in `molspotter`. They are respectively trained on:
- [ChEMBL](https://doi.org/10.1093/nar/gky1075) (version 32)
- [ExCAPE-DB](https://doi.org/10.1186/s13321-017-0203-5) [(Zenodo)](https://doi.org/10.5281/zenodo.675987)
- [Papyrus](https://doi.org/10.1186/s13321-022-00672-x) [(version 05.6)](https://doi.org/10.5281/zenodo.7373213)

### API
One can easily load a spotter as follows:

```python
from molspotter import SillyMolSpotter

sms = SillyMolSpotter.from_pretrained('chembl')
```

New molecules can easily be scored using this spotter instance:

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('C=C=C1CN(C)CCN1Cc5ccc(C(=O)Nc4ccc(C)c(Nc3nccc(c2cccnc2)n3)c4)cc5')
sms.score_mol(mol)
# 1
```
The binary silly bit score is defined as 1 if the molecule contains any substructure not observed in the reference set and 0 otherwise.
The continuous score is defined as $score = \frac{n_{silly bits}}{n_{on bits}}$.
A value of 0.0 means that the molecule does not contain any silly bit, while a value of 1.0 means the molecule is only made of silly bits.


Silly bits can also be displayed:
```python
img = sms.show_mol(mol)
```

![Silly bits highlighted in the molecular structure](images/silly_bits.png)


One can create their own spotter using the following:
 ```python
# Instantiate a spotter with radius of Morgan fingerprint value of 2
sms = SillyMolSpotter(fp_radius=2)

sms.add_file('PATH_TO_MOLECULAR_FILE')
```

Finally, a spotter can be saved and loaded back from disk:

```python
sms.save('PATH_TO_SAVE_SPOTTER_TO')

sms2 = SillyMolSpotter.from_file('PATH_TO_SAVED_SPOTTER')
```

### Command-line interface

- Create a spotter
```commandline
molspotter create -o examples/chembl_drugs.sp -i examples/chembl_drugs.sd
```

- Score with a spotter
```commandline
molspotter score -s examples/chembl_drugs.sp -i examples/chembl_drugs.sd
molspotter score -s chembl -i examples/chembl_drugs.sd
```
