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
    parser.add_argument("--metadata1", required=True, help="Metadata file from NextStrain")
    parser.add_argument("--metadata2", required=False, help="Custom lab metadata file")
    parser.add_argument("--output1", required=True, help="Filtered metadata file")
    parser.add_argument("--output2", required=True, help="TSV file for renaming virus IDs")
    parser.add_argument("--output3", required=True, help="Reformatted, final FASTA file")
    args = parser.parse_args()

    genomes = args.genomes
    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output1 = args.output1
    output2 = args.output2
    output3 = args.output3

    # genomes = path + 'temp_sequences.fasta'
    # metadata1 = path + 'metadata_nextstrain.tsv'
    # metadata2 = path + 'COVID-19_sequencing.xlsx'
    # output1 = path + 'metadata_filtered.tsv'
    # output2 = path + 'rename.tsv'
    # output3 = path + 'sequences.fasta'

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
    dfN = pd.read_csv(metadata1, encoding='utf-8', sep='\t', dtype='str')
    try:
        dfN = dfN[['strain', 'gisaid_epi_isl', 'genbank_accession', 'date', 'country', 'division', 'location',
                      'region_exposure', 'country_exposure', 'division_exposure', 'originating_lab', 'submitting_lab',
                   'authors']]
        dfN.insert(4, 'iso', '')
    except:
        pass
    dfN['update'] = ''
    dfN.fillna('', inplace=True)
    lColumns = dfN.columns.values  # list of column in the original metadata file

    # Lab genomes metadata
    dfL = pd.read_excel(metadata2, index_col=None, header=0, sheet_name='Amplicon_Sequencing',
                        # 'sheet_name' must be changed to match the Excel sheet name
                        converters={'Sample-ID': str, 'Collection-date': str,
                                    'Update': str})  # this need to be tailored to your lab's naming system
    dfL.fillna('', inplace=True)
    dfL.set_index("Sample-ID", inplace=True)


    dHeaders = {}
    notFound = []
    lstNewMetadata = []
    found = []
    lab_label = {}
    for id in sequences.keys():
        # check nextstrain metadata first
        dRow = {}
        if id in dfN['strain'].to_list() and 'Yale-' not in id:
            fields = {column: '' for column in lColumns}
            row = dfN.loc[lambda dfN: dfN['strain'] == id]

            strain = row.strain.values[0]
            country = row.country.values[0]
            division = row.division.values[0]
            location = row.location.values[0]

            country_exposure = row.country_exposure.values[0].strip()
            if country != country_exposure:  # ignore travel cases
                continue

            division_exposure = row.division_exposure.values[0].strip()
            if division != division_exposure:  # ignore travel cases
                continue

            if len(country) < 2:
                row.country.values[0] = ''
            if len(division) < 2:
                row.division.values[0] = country

            iso = get_iso(country)
            row.iso.values[0] = iso  # needed for exporting a renaming file
            date = row.date.values[0]
            header = '|'.join([strain, iso, division.replace(' ', '-'), date])
            dHeaders[strain] = header

            lValues = row.values[0]
            for field, value in zip(fields.keys(), lValues):
                # print(id)
                if value in ['', np.nan, None]:
                    value = ''

                fields[field] = value
            if country == '':
                continue

            dRow[id] = fields
            found.append(strain)
            print('Exporting metadata for ' + id)

        # check lab's metadata otherwise
        if id not in dRow.keys():
            # check lab metadata
            if 'Yale-' in id:
                try:
                    id = id.split('/')[1][3:]
                    lab_label[id] = ''
                except:
                    if id not in lab_label.keys():
                        lab_label[id] = ''

                if id in dfL.index:
                    fields = {column: '' for column in lColumns}
                    row = dfL.loc[id]
                    if row['State'] == '':
                        code = 'CT'  # change this line to match the acronym of the most likely state of origin if the 'State' field is unknown
                    else:
                        code = row['State']
                    strain = 'USA/' + code + '-' + id + '/2020'  # change this line to match the country of origin (alpha-3 ISO code)

                    if strain not in found:
                        gisaid_epi_isl = ''
                        genbank_accession = ''
                        if len(str(row['Collection-date'])) > 1:
                            date = row['Collection-date'].split(' ')[0].replace('.', '-').replace('/', '-')
                        else:
                            date = ''

                        country = row['Country']

                        division = row['Division']
                        if row['Division'] in ['', '?']:
                            code = 'Connecticut'  # change this line to match the most likely state of origin if the 'Division' field is unknown
                        else:
                            code = row['Division']

                        if row['Location'] in ['', '?', 'N/A']:
                            location = ''
                        else:
                            location = str(row['Location'])

                        region_exposure = ''
                        country_exposure = ''
                        iso = 'USA'  # change this line to match the country of origin (alpha-3 ISO code)
                        division_exposure = ''
                        try:
                            length = str(len(sequences[strain]))
                        except:
                            length = str(len(sequences[id]))
                        host = row['Host']
                        originating_lab = row['Source']
                        submitting_lab = 'Grubaugh Lab - Yale School of Public Health'  # change this line to match you lab's name
                        authors = 'Fauver et al'  # change this line to match you lab's main author's name
                        update = 'Update' + str('0' * (2 - len(row['Update']))) + row['Update']

                        lValues = [strain, gisaid_epi_isl, genbank_accession, date, iso, country, division,
                                   location, region_exposure, country_exposure, division_exposure,
                                   originating_lab, submitting_lab, authors, update]

                        header = '|'.join([strain, country, division.replace(' ', '-'), date])
                        dHeaders[strain] = header

                        for field, value in zip(fields.keys(), lValues):
                            fields[field] = value
                        dRow[id] = fields
                        found.append(strain)
                        if id in lab_label.keys():
                            lab_label[id] = strain
                    else:
                        continue

            else:  # Assign 'NA' if no metadata is available
                header = '|'.join([id, 'NA', 'NA', 'NA', 'NA'])
                dHeaders[id] = header
                notFound.append(id)
        lstNewMetadata = lstNewMetadata + list(dRow.values())


    # write new metadata files
    outputDF = pd.DataFrame(lstNewMetadata, columns=list(lColumns))
    outputDF.to_csv(output1, sep='\t', index=False)

    # write renaming file
    with open(output2, 'w') as outfile2:
        # export new metadata lines
        for id, header in dHeaders.items():
            outfile2.write(id + '\t' + header + '\n')
        for id in notFound:
            print('\t* Warning! No metadata found for ' + id)

        if len(notFound) > 0:
            print('\nPlease check for inconsistencies (see above).')

    # write sequence file
    exported = []
    with open(output3, 'w') as outfile3:
        # export new metadata lines
        for id, sequence in sequences.items():
            if 'Yale' in id:
                if lab_label[id] not in exported:
                    entry = '>' + lab_label[id] + '\n' + sequence + '\n'
                    outfile3.write(entry)
                    print('* Exporting newly sequenced genome and metadata for ' + id)
                    exported.append(lab_label[id])
            else:
                if id not in exported:
                    entry = '>' + id + '\n' + sequence + '\n'
                    outfile3.write(entry)
                    exported.append(id)

print('\nMetadata file successfully reformatted and exported!\n')
