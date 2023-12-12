"""
python [input dir] [metadata csv] [output csv]

input dir: directory with star-fusion directories for each sample
star-fusion directories should have star-fusion.fusion_predictions.abridged.coding_effect.tsv
file

"""

import os
import sys
import pandas as pd

def concat_fusions(files):
    fusions = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        accession = os.path.basename(os.path.dirname(f))
        sample_id = accession.split('-')[0]
        print(f"Getting fusions for {sample_id}")
        df['sample_id'] = sample_id
        fusions.append(df)
    return pd.concat(fusions)

def main():
    star_dirs_fp = sys.argv[1]
    metadatafp = sys.argv[2]
    outfp = sys.argv[3]

    fusion_fn = 'star-fusion.fusion_predictions.abridged.coding_effect.tsv' 

    sample_dirs = os.listdir(star_dirs_fp)
    fusion_files = [os.path.join(star_dirs_fp, x, fusion_fn) for x in sample_dirs]
    # fusion_files = fusion_files[:3]

    metadata = pd.read_csv(metadatafp)
    metadata['sample_id'] = metadata['patient']
    metadata = metadata[['sample_id', 'condition']]
    
    fusions = concat_fusions(fusion_files)
    # fusions.to_csv('fusions.csv', index=False)
    fusions = fusions.merge(metadata)
    fusions.to_csv(outfp, index=False)

if __name__ == '__main__':
    main()
