dc = {}
for i in df.iterrows():
    ids = i[1].Protein_id.split(';')
    for j in ids:
        k = j.split('-')[0]
        if len(df2[df2.Protein_id == k]) > 0:
            dc[i[0]] = list(df2[df2.Protein_id == k].Protein_id)[0]
            break
