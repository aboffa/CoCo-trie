K = 31
with open('dna', 'r') as file:
    data = file.read().replace('\n', '')
    res = set(data[i:i+K] for i in range(len(data)-K+1))
    res = sorted(list(res))

    with open('dna-'+str(K)+'-mer.txt', mode='w') as f:
        f.write('\n'.join(res))
