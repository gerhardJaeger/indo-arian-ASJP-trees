f = open('language_coding.txt','r')

text = f.read()

f.close()

text = text.split('\n\n\n')

states = [t.split('\n') for t in text]
states = sorted(states)

f = open('states.csv','w')

for t in states:
    print(t[0].replace(' ','_')+'<'+t[1]+'>',t[3],file=f)


f.close()

f = open('states_appendix.tex','w')

for l in states:
    print('\paragraph{'+l[0]+' ['+l[1]+'] '+l[3]+'}',file=f)
    print(l[2],file=f)


f.close()
