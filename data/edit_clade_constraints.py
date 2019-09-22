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