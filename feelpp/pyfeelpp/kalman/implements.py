def balancedpartition(nb_data,nb_procs):
    part = np.zeros(nb_procs, dtype=np.int8)
    for i in range(nb_procs):
        part[i] = int(round(nb_data/(nb_procs-i)))
        nb_data -= part[i]
    return part

def displacements(partition):
    len = np.size(partition)
    disp = np.zeros(len)
    for i in range(1,len):
        disp[i] = disp[i-1] + partition[i-1]
    return disp
    
