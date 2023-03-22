# -*- coding: utf-8 -*-

import os

import click

from . import SillyMolSpotter, MolSupplier


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    pass


@main.command(help='Create a spotter.', context_settings=CONTEXT_SETTINGS)
@click.option('-o', '--outfile', type=str, required=True,
              default=None, nargs=1, show_default=True, metavar='OUTFILE',
              help='Path to save the spotter to.')
@click.option('-r', '--radius', required=False, default=1, nargs=1,
              metavar='RADIUS', help='Radius of the Morgan fingerprinter to be used to identify silly bits.')
@click.option('-i', '--input', required=True, multiple=True, show_default=True, metavar='INFILE',
              help='Input file(s) to be parsed and added to the spotter.')
def create(outfile, radius, input):
    if isinstance(input, tuple):
        input = list(input)
    sms = SillyMolSpotter(radius)
    for file in input:
        sms.add_file(file)
    sms.save(outfile)


@main.command(help='Scores molecules with a spotter.', context_settings=CONTEXT_SETTINGS)
@click.option('-s', '--spotter', type=str, required=True,
              default=None, nargs=1, show_default=True, metavar='SPOTTER',
              help='Path to the saved the spotter or one of \'ChEMBL\', \'ExCAPE\' or \'Papyrus\'.')
@click.option('-i', '--input', required=True, multiple=True,
              metavar='INPUT', help='Input file(s) containing molecules to be scored.')
@click.option('-o', '--output', required=False, nargs=1,
              metavar='OUTPUT', help='Output file in which scores will be stored.')
def score(spotter, input, output):
    if output is not None and not os.access(os.path.dirname(output), os.W_OK):
        raise IOError(f'cannot write file to {output}')
    # try:
    sms = SillyMolSpotter.from_pretrained(spotter)
    # except:
    #     sms = SillyMolSpotter.from_file(spotter)
    if isinstance(input, tuple):
        input = list(input)
    if output is not None:
        output_handle = open(output, 'w')
    for file in input:
        with MolSupplier(file) as supplier:
            for mol in supplier:
                score = sms.score_mol(mol)
                if output is not None:
                    output_handle.write(str(score) + '\n')
                else:
                    print(score)
    if output is not None:
        output_handle.close()
