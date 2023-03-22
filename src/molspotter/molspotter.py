# -*- coding: utf-8 -*-

from __future__ import annotations
import os
import json
import lzma
from typing import List, Union

from rdkit import Chem
from rdkit.rdBase import BlockLogs
from rdkit.Chem import AllChem, Draw
from PIL.ImageFile import ImageFile

from .io import MolSupplier


CHEMBL_SPOTTER = os.path.join(os.path.dirname(__file__), 'spotters', 'chembl_32.sp')
EXCAPE_SPOTTER = os.path.join(os.path.dirname(__file__), 'spotters', 'excapedb.sp')
PAPYRUS_SPOTTER = os.path.join(os.path.dirname(__file__), 'spotters', 'papyrus_05.6.sp')


class SillyMolSpotter:
    def __init__(self, fp_radius):
        """Initialize a spotter to identify silly molecules.

        :param fp_radius: radius of the Morgan fingerprinter
        """
        self.count_dict = {}
        self.radius = fp_radius

    @classmethod
    def from_file(cls, filepath: str) -> SillyMolSpotter:
        """Load a saved spotter."""
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f'file does not exist {filepath}')
        with lzma.open(filepath) as infile:
            data = json.load(infile)
        sms = SillyMolSpotter(data['spotter_radius'])
        sms.count_dict = {int(x): int(y) for x, y in data['spotter_bits'].items()}
        return sms

    @classmethod
    def from_pretrained(cls, name) -> SillyMolSpotter:
        """Load a pretrained spotter.

        :param name: one of 'ChEMBL', 'ExCAPE' or 'Papyrus'
        """
        mapping = {'chembl': CHEMBL_SPOTTER,
                   'excape': EXCAPE_SPOTTER,
                   'papyrus': PAPYRUS_SPOTTER}
        if name.lower() not in mapping.keys():
            raise ValueError('Name of a pretrained spotter may only be one of \'ChEMBL\', \'ExCAPE\' or \'Papyrus\'')
        sms = SillyMolSpotter.from_file(mapping[name])
        return sms

    def add_mol(self, mol: Union[Chem.Mol, List[Chem.Mol]]) -> None:
        """Add sound molecule(s) to the reference set of the spotter.

        :param mol: a RDKit molecule or list of molecules
        """
        if not isinstance(mol, list):
            mol = [mol]
        block = BlockLogs()
        for mol_ in mol:
            if mol_ is not None:
                try:
                    fp = AllChem.GetMorganFingerprint(mol, self.radius)
                    for k, v in fp.GetNonzeroElements().items():
                        self.count_dict[k] = self.count_dict.get(k, 0) + v
                except:
                    pass
        del block

    def add_smiles(self, smiles: Union[str, List[str]]) -> None:
        """Add SMILES of molecule(s) to the reference set of the spotter.

        :param smiles: SMILES of a molecules or of a list of molecules.
        """
        if not isinstance(smiles, list):
            smiles = [smiles]
        block = BlockLogs()
        for notation in smiles:
            try:
                mol = Chem.MolFromSmiles(notation)
                if mol is not None:
                    fp = AllChem.GetMorganFingerprint(mol, self.radius)
                    for k, v in fp.GetNonzeroElements().items():
                        self.count_dict[k] = self.count_dict.get(k, 0) + v
            except:
                pass
        del block

    def add_file(self, filepath: str) -> None:
        """Add the molecular content of a file to the reference set of the spotter.

        :param filepath: path to a file containing molecular structures.
        """
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f'file cannot be found: {filepath}')
        file_size = bool(os.path.getsize(filepath) // 2 ** 20)  # file is more than 1MB
        block = BlockLogs()
        with MolSupplier(filepath, show_progress=file_size) as supplier:
            for mol in supplier:
                if mol:
                    try:
                        fp = AllChem.GetMorganFingerprint(mol, self.radius)
                        for k, v in fp.GetNonzeroElements().items():
                            self.count_dict[k] = self.count_dict.get(k, 0) + v
                    except:
                        pass
        del block

    def score_mol(self, mol: Union[Chem.Mol, List[Chem.Mol]]) -> Union[float, List[float]]:
        """Score how silly the input molecule(s) is/are.

        :param mol: RDKit molecule or list of molecules
        """
        if not len(self.count_dict):
            raise ValueError('this spotter does not contain any bit, please add molecules before scoring')
        if not isinstance(mol, list):
            mol = [mol]
        scores = []
        for mol_ in mol:
            if mol_ is not None:
                fp = AllChem.GetMorganFingerprint(mol_, self.radius)
                on_bits = fp.GetNonzeroElements().keys()
                silly_bits = [x for x in [self.count_dict.get(x)
                                          for x in on_bits]
                              if x is None]
                score = len(silly_bits) / len(on_bits)
            else:
                score = 1
            scores.append(score)
        if len(scores) == 1:
            return scores[0]
        return scores

    def save(self, filepath: str) -> None:
        """Save a spotter to disk.

        :param filepath: path on disk to save the spotter to.
        """
        if not os.access(os.path.dirname(filepath), os.W_OK):
            raise IOError(f'cannot write file to {filepath}')
        with lzma.open(filepath, mode='wt', preset=9 | lzma.PRESET_EXTREME) as outfile:
            json.dump({'spotter_radius': self.radius,
                       'spotter_bits': self.count_dict},
                      outfile)

    def show_mol(self, mol: Chem.Mol, raise_error: bool = False, molsPerRow=3,
             subImgSize=(500, 300), useSVG=False, returnPNG=False, **kwargs) -> Union[str, ImageFile]:
        block = BlockLogs()
        if mol:
            # Extract bits information
            bi = {}
            _ = AllChem.GetMorganFingerprint(mol, radius=self.radius, bitInfo=bi, useCounts=False)
            # Collect atom and bond environments for drawing
            drawings = []
            for hash, atom_abb_envs in bi.items():
                skip = False
                atom_env = set()
                bond_env = []
                for atom, radius in atom_abb_envs:
                    if self.count_dict.get(hash) is not None:
                        skip = True
                        continue
                    bond_env.extend(Chem.FindAtomEnvironmentOfRadiusN(mol, self.radius, atom))
                    for bond in bond_env:
                        atom_env.add(mol.GetBondWithIdx(bond).GetBeginAtomIdx())
                        atom_env.add(mol.GetBondWithIdx(bond).GetEndAtomIdx())
                if not skip:
                    legend = f'{hash}'
                    drawings.append((tuple(atom_env), tuple(bond_env), legend))
            # Silly bits found
            if len(drawings):
                atom_highlights, bond_highlights, legends = zip(*drawings)
                return Draw.MolsToGridImage([mol] * len(drawings),
                                            molsPerRow=molsPerRow,
                                            subImgSize=subImgSize,
                                            legends=legends,
                                            highlightAtomLists=atom_highlights,
                                            highlightBondLists=bond_highlights,
                                            useSVG=useSVG,
                                            returnPNG=returnPNG,
                                            **kwargs)
            return 'No silly bit'
            # return Draw.MolToImage(mol, highlightAtoms=atoms_envs, highlightBonds=bonds_envs, legend='3026373406')
        elif raise_error:
            del block
            raise TypeError('smiles could not be parsed')
        del block
