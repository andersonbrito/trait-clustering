# -*- coding: utf-8 -*-

import pycountry_convert as pyCountry
import pycountry
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genomes", required=True, help="FASTA file genomes to be used")
    parser.add_argument("--metadata", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata = args.metadata
    output1 = args.output1
    output2 = args.output2


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_immune/nextstrain/run1_test/pre-analyses/'
    # genomes = path + 'temp_sequences.fasta'
    # metadata = path + 'metadata_nextstrain.tsv'
    # output1 = path + 'metadata_filtered.tsv'
    # output2 = path + 'sequences.fasta'


    print('\n### Parsing fasta file\n')
    # create a dict of existing sequences
    sequences = {}
    for fasta in SeqIO.parse(open(genomes), 'fasta'):  # as fasta:
        id, seq = fasta.description, fasta.seq
        if id not in sequences.keys():
            sequences[id] = str(seq)

    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # nextstrain metadata
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype='str')
    dfN.insert(4, 'iso', '')
    dfN.fillna('', inplace=True)

    selected = list(sequences.keys())
    dfN = dfN[dfN['strain'].isin(selected)]
    dfN = dfN[~dfN['country'].isin(['', 'unknown', 'na'])]

    # fix exposure
    for exposure_column in ['region_exposure', 'country_exposure', 'division_exposure']:
        for idx, row in dfN.iterrows():
            level = exposure_column.split('_')[0]
            if dfN.loc[idx, exposure_column].lower() in ['', 'unknown', 'na']:
                dfN.loc[idx, exposure_column] = dfN.loc[idx, level]

    dfN['iso'] = dfN['country_exposure'].map(get_iso)
    notFound = [id for id in selected if id not in dfN['strain'].to_list()]

    print('### Exporting metadata\n')
    # write new metadata files
    dfN.to_csv(output1, sep='\t', index=False)

    # write sequence file
    exported = []
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, sequence in sequences.items():
            if id not in exported and id not in notFound:
                entry = '>' + id + '\n' + sequence + '\n'
                outfile2.write(entry)
                exported.append(id)

    c = 1
    if len(notFound) > 0:
        print('\tNo metadata found for these entries:\n')
        for entry in notFound:
            print('\t' + str(c) + '. ' + entry)
            c += 1

    print('\n\nSequence and Metadata files successfully exported!\n')
    print('\t- ' + str(len(notFound)) + ' genomes have no metadata, and were removed.')
    print('\t- ' + str(len(exported)) + ' genomes were correctly exported.\n')
