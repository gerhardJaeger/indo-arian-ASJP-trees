from ete3 import Tree


ctree = Tree('IIr_glotto_tree.newick',format=1)


#nodes_to_keep = [n.name for n in ctree.get_descendants() if not n.is_leaf()]
nodes_to_keep = ['PIA', 
                 'PIr', 
                 #'PNonSkt', 
                 #'PnonAv', 
                 'PNonDard', 
                 #'PDard', 
                 'wakh|khot', 
                 #'PnonSaka', 
                 #'PContIA', 
                 'dhiv|sinh', 
                 #'PnonShin', 
                 'sari|shug|ishk', 
                 #'ormu|para|pash', 
                 #'PnonPash', 
                 'panj|sera|sind', 
                 #'PNucIA', 
                 #'PWIA', 
                 #'khow|kala', 
                 #'phal|torw|indu', 
                 'shug|ishk', 
                 'ormu|para', 
                 #'PEI1', 
                 #'PWIr', 
                 'panj|sera', 
                 'PEIA', 
                 'PCentWIA', 
                 'guja|saur', 
                 'mara|konk', 
                 #'torw|indu', 
                 #'khwa|sogd|yagn', 
                 #'west|cent|zaza', 
                 'PSWIr', 
                 #'Caspian', 
                 #'assa|beng|oriy', 
                 #'PCentEIA', 
                 'CentIA', 
                 'PRaj', 
                 'sogd|yagn', 
                 #'cent|zaza', 
                 'OP', 
                 #'PLuri', 
                 #'gila|sang|esht', 
                 #'beng|oriy', 
                 #'awad|bhoj|mait', 
                 'vlax|doma|doma', 
                 #'marw|dhun', 
                 #'mewa|bagr', 
                 'MP', 
                 #'jude|luri|bakh', 
                 #'gila|sang', 
                 #'awad|bhoj', 
                 #'vlax|doma', 
                 #'midd|dari', 
                 #'luri|bakh'
                 ]


for n in ctree.get_descendants():
    if not n.is_leaf():
        if n.name not in nodes_to_keep:
            #print (n.name)
            n.delete()



ctree.write(format=9,outfile='constraints2.tre')


clade_calib_anc = {
    'PIA':('sans1269','Unif(3400,3000)'),
    'OP':('oldp1254','Unif(2500,2300)'),
    'MP':('pahl1241','Unif(1800,1400)')
}


clade_calib = {
    'PIIr':'Unif(3900,3600)', #post BMAC breakup
    'PIr':'Unif(3500,3300)'   #Yaz I
}


f = open('node_calibration.txt','w')

print("""#each block gives (1) the clade name
#(2) tip which should have a very short branch, if relevant
#(3) the distribution over the age of the clade
#(4) the tips in the clade
""",file=f)

for k in clade_calib_anc.keys():
    print(k,file=f)
    print(clade_calib_anc[k][0],file=f)
    print(clade_calib_anc[k][1],file=f)
    print(' '.join(ctree.search_nodes(name=k)[0].get_leaf_names()),file=f)
    print('\n',file=f)


for k in clade_calib.keys():
    print(k,file=f)
    print(clade_calib[k],file=f)
    print(' '.join(ctree.search_nodes(name=k)[0].get_leaf_names()),file=f)
    print('\n',file=f)


f.close()