import pandas as pd
import os
import numpy as np
from ete3 import Tree

def nexCharOutput(chMtx, outfile, datatype='STANDARD'):
    names = chMtx.index
    with open(outfile, 'w') as f:
        f.write('#NEXUS\n\n')
        f.write('BEGIN DATA;\n')
        f.write('DIMENSIONS ntax=' + str(len(chMtx)) +
                ' NCHAR='+str(len(chMtx.T))+';\n')
        f.write('FORMAT DATATYPE=' + datatype +
                ' GAP=? MISSING=- interleave=yes;\n')
        f.write('MATRIX\n\n')
        txLgth = max(map(len, names))
        for i in range(len(chMtx)):
            f.write(names[i].ljust(txLgth+2))
            for ch in chMtx.values[i]:
                if ch == -1:
                    ch = '-'
                else:
                    ch = str(ch)
                f.write(ch)
            f.write('\n')
        f.write('\n;\n\nEND;\n')
        f.close()



asjpCC = pd.read_csv('../data/asjp18Clustered.csv',
                     na_filter=False)
asjp = pd.read_csv('../data/dataset.tab', sep='\t',
                   na_filter=False)

sounds = np.unique(list(''.join(asjpCC.simplified.values)))

concepts = asjpCC.concept.unique()

asjpCC.doculect = [x.replace('-', '_').replace("'", "")
                   for x in asjpCC.doculect.values]
asjp['names'] = [x.replace('-', '_').replace("'", "")
                 for x in asjp.names.values]

asjp['longname'] = ['.'.join(x).replace('-', '_').replace("'", "")
                    for x in asjp[['wls_fam', 'wls_gen', 'names']].values]


asjp = asjp[~asjp.iso.isnull()]


asjp['nEntries'] = asjp[concepts].apply(lambda x:
                                        sum(~x.isnull()), axis=1)

asjp = asjp[asjp.names.isin(asjpCC.doculect)]

asjp.sort_values('nEntries',
                 ascending=False,
                 inplace=True)

metadata = pd.read_csv('../data/languages.csv')


metadata['doculect'] = [x.replace('-','_') for x in metadata.ID.values]


indoiranian = pd.read_csv('../data/indo_iranian_asjp_taxa.tsv',
                          sep='\s+', names=['glottocode', 'doculect', 'comment'],
                          header=None)


def longname2glottocode(l):
    if l not in asjp.longname.values:
        return np.nan
    d = asjp[asjp.longname == l].names.values[0]
    if d not in metadata.doculect.values:
        return np.nan
    return metadata[metadata.doculect == d].Glottocode.values[0]

l2g = pd.Series({l: longname2glottocode(l)
                 for l in asjp.longname.values})


asjp['glottocode'] = l2g[asjp.longname.values].values

asjp = asjp.drop_duplicates('glottocode')

def doculect2glottocode(d):
    if d not in metadata.doculect.values:
        return np.nan
    return metadata[metadata.doculect == d].Glottocode.values[0]

d2g = pd.Series({l: doculect2glottocode(l)
                 for l in asjp.names.values})

asjpCC['glottocode'] = d2g[asjpCC.doculect.values].values

asjpCC = asjpCC[asjpCC.doculect.isin(asjp.names)]

asjpCC = asjpCC[asjpCC.glottocode.isin(indoiranian.glottocode.values)]

ctree = Tree('../data/constraint3.tre')
ctree.ladderize()

taxa = indoiranian.glottocode.unique()

ccMtx = pd.DataFrame(index=taxa)
for c in concepts:
    cData = asjpCC[asjpCC.concept == c].copy()
    cMtx = pd.crosstab(cData.glottocode, cData.cClass)
    cMtx[cMtx > 1] = 1
    cMtx = cMtx.reindex(taxa, fill_value='-')
    ccMtx = pd.concat([ccMtx, cMtx], axis=1)

scMtx = pd.DataFrame(index=taxa)
for c in concepts:
    cData = asjpCC[asjpCC.concept == c]
    cTaxa = cData.glottocode.unique()
    cWords = pd.Series([''.join(cData[cData.glottocode ==
                                      l].simplified.values)
                        for l in cTaxa],
                       index=cTaxa)
    cMtx = pd.DataFrame([[int(s in cWords[l]) for s in sounds]
                         for l in cTaxa],
                        index=cTaxa,
                        columns=[c+':'+s
                                 for s in sounds]).reindex(taxa,
                                                           fill_value='-')
    scMtx = pd.concat([scMtx, cMtx], axis=1)

n, m = ccMtx.shape[1], scMtx.shape[1]
cc_sc = pd.concat([ccMtx, scMtx], axis=1)

nexCharOutput(cc_sc, 'indoiranian.nex', datatype='restriction')

nodes = np.array([nd for nd in ctree.get_descendants()
                  if not nd.is_root() and not nd.is_leaf()])

for i, nd in enumerate(nodes):
    nd.name = 'clade'+str(i+1).rjust(2, '0')

def nname(x):
    if x.is_leaf():
        return '"'+x.name+'"'
    else:
        return x.name


rev = ""
for nd in reversed(nodes):
    rev += nd.name + " = clade("
    rev += ', '.join([nname(x) for x in nd.get_children()])
    rev += ")\n"
rev += 'constraints = ['+', '.join([nname(nd) for nd in reversed(nodes)])
rev += ']\n'

with open('constraints.Rev', 'w') as f:
    f.write(rev)
