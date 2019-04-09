###This performs the topology profiling on a set of introgressed blocks that are annotated
###with mutation motifs.

###Mutation motifs can be retrieved using vcftools or bcftools by e.g. cutting down VCFs
###to the 4 relevant individuals [Target,Deni,Nean,HumanOutgroup] where HumanOutgroup is
###a human sequence (for our introgression analysis I used a Russian sample,
###LP6005441-DNA_G10); trimming down to the introgressed regions; and counting the
###occurence of each of the 16 possible mutation motifs.

###The function calculateMotifs_BEDonBCF.py can be used to retrieve input files.

import numpy as np

def topology_analysis(topology_files_in = ['./example_topoMotifCounts_CDNH.txt']):
    """
    This is an extended topology analysis.
    It reads in a set of full-topology-annotated chunks. It converts to mismatches, and calculates pairwise distances between each chunk set.
    It then tries to assign consistency of each chunk with the 15 possible coalescent trees.
    """
    #The format is [name, motifs, mismatch, pairwise, topology support]
    #I assume we are [Chunk,Deni,Nean,Human]
    #(Note that this motif order is different from the published analysis)
    #          0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
    #motifs:   0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
    #          0    1    2    3    4    5    6    7    8    9
    #pairwise: CDen CNea CHum CAnc DNea DHum DAnc NHum NAnc HAnc
    #topology: (H,(N,(D,C))), (H,(D,(N,C))), (H,(C,(D,N))),
    #          (N,(H,(D,C))), (N,(D,(H,C))), (N,(C,(D,H))),
    #          (D,(N,(H,C))), (D,(H,(N,C))), (D,(C,(H,N))),
    #          (C,(N,(D,H))), (C,(D,(N,H))), (C,(H,(D,N))),
    #          ((H,N),(D,C)), ((H,D),(N,C)), ((H,C),(D,N))
    topo_order = ['(H,(N,(D,C)))', '(H,(D,(N,C)))', '(H,(C,(D,N)))', '(N,(H,(D,C)))', '(N,(D,(H,C)))', '(N,(C,(D,H)))', '(D,(N,(H,C)))', '(D,(H,(N,C)))', '(D,(C,(H,N)))', '(C,(N,(D,H)))', '(C,(D,(N,H)))', '(C,(H,(D,N)))', '((H,N),(D,C))', '((H,D),(N,C))', '((H,C),(D,N))']

    inconsistency_threshold = 0.10
    branch_length_equality_factor = 1.5 #If 1.5 (or 2/3.), only reject if a pair of equal branches are have a 1:1.5 or 1.5:1 ratio
    min_unmasked = 20000
    topology_focus = [[0], [12]] #Optional assessment of mismatch in specific subset of topologies. Here I check divergence in topologies (H,(N,(D,C))) [0] and ((H,N),(D,C)) [12].
    #topology_focus = []

    set_names = ['example_topos']
    ind_column = 4
    chromosomes_to_exclude = []
    
    mismatch_keys = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']
    pairwise_keys = [[0,1],[0,2],[0,3],[0,'anc'],[1,2],[1,3],[1,'anc'],[2,3],[2,'anc'],[3,'anc']]
    #Each topology has four consistency tests:

    #1. That the mismatches are in the correct order. Here, the notation is [0,...,9] with [[a],[b],[c]] indicating a < b < c and [[a,b],[c]] indicating \bar{a,b} < x
    #2. The the (average) total branch lengths are too. Here, the notation is [0,...,15] with [[[a]],[[b,d],[c,d]]] indicating a < \bar{b+d,c+d}. NB that the outer list is averaging and the inner list is summation.
    #3. That the number of inconsistent motifs is low. Here, the notation [0,...,15] with the list indicating the motifs (see mismatch_keys) that are inconsistent with the proposed topology. We require that under inconsistency_threshold motifs are inconsistent with a topology to accept that topology.
    #4. That branches that should be of equal length given the topology are approximately equal. The notation is [0,...,15] indicating motifs (see mismatch_keys) that contribute to different branches under the assumption of the proposed topology. An input of [[[a],[b]], [[c,d],[e]]] would have two conditions - that abs(ln(a/b)) < abs(ln(branch_length_equality_factor)) and that abs(ln((c+d)/e)) < abs(ln(branch_length_equality_factor)). I.e. [[[motifs contributing to branch 1],[motifs contributing to branch 2 that should be equal to branch 1]], [[motifs contributing to branch 3], [motifs contributing to branch 4, that should be equal to branch 3]], ...other conditions...]
    
    topo_conditions = [[[[[0],[4],[2]],[[0],[4],[5]],[[0],[4],[7]],[[0],[1],[2]],[[0],[1],[5]],[[0],[1],[7]]], [[[[4],[8]],[[2]]], [[[2],[4,12],[8,12]],[[1]]]], [3,5,7,9,11,13,6,10],   [[[4],[8]], [[4,12],[2]], [[8,12],[2]], [[4,12,14],[1]], [[8,12,14],[1]]]],
                       [[[[1],[0],[2]],[[1],[0],[5]],[[1],[0],[7]],[[1],[4],[2]],[[1],[4],[5]],[[1],[4],[7]]], [[[[2],[8]],[[4]]], [[[4],[2,10],[8,10]],[[1]]]], [3,5,7,9,11,13,6,12],   [[[2],[8]], [[2,10],[4]], [[8,10],[4]], [[2,10,14],[1]], [[8,10,14],[1]]]],
                       [[[[4],[0],[2]],[[4],[0],[5]],[[4],[0],[7]],[[4],[1],[2]],[[4],[1],[5]],[[4],[1],[7]]], [[[[2],[4]],[[8]]], [[[8],[2,6],[4,6]],[[1]]]],   [3,5,7,9,11,13,10,12],  [[[2],[4]], [[2,6],[8]], [[4,6],[8]], [[2,6,14],[1]], [[4,6,14],[1]]]],
                       
                       [[[[0],[2],[1]],[[0],[2],[4]],[[0],[2],[7]],[[0],[5],[1]],[[0],[5],[4]],[[0],[5],[7]]], [[[[4],[8]],[[1]]], [[[1],[4,12],[8,12]],[[2]]]], [3,6,7,10,11,14,5,9],   [[[4],[8]], [[4,12],[1]], [[8,12],[1]], [[4,12,13],[2]], [[8,12,13],[2]]]],
                       [[[[2],[0],[1]],[[2],[0],[4]],[[2],[0],[7]],[[2],[5],[1]],[[2],[5],[4]],[[2],[5],[7]]], [[[[1],[8]],[[4]]], [[[4],[1,9],[8,9]],[[2]]]],   [3,6,7,10,11,14,5,12],  [[[1],[8]], [[1,9],[4]], [[8,9],[4]], [[1,9,13],[2]], [[8,9,13],[2]]]],   
                       [[[[5],[0],[1]],[[5],[0],[4]],[[5],[0],[7]],[[5],[2],[1]],[[5],[2],[4]],[[5],[2],[7]]], [[[[1],[4]],[[8]]], [[[8],[1,5],[4,5]],[[2]]]],   [3,6,7,10,11,14,9,12],  [[[1],[4]], [[1,5],[8]], [[4,5],[8]], [[1,5,13],[2]], [[4,5,13],[2]]]],
                       
                       [[[[2],[1],[0]],[[2],[1],[4]],[[2],[1],[5]],[[2],[7],[0]],[[2],[7],[4]],[[2],[7],[5]]], [[[[1],[8]],[[2]]], [[[2],[1,9],[8,9]],[[4]]]],   [5,6,7,12,13,14,3,10],  [[[1],[8]], [[1,9],[2]], [[8,9],[2]], [[1,9,11],[4]], [[8,9,11],[4]]]],
                       [[[[1],[2],[0]],[[1],[2],[4]],[[1],[2],[5]],[[1],[7],[0]],[[1],[7],[4]],[[1],[7],[5]]], [[[[2],[8]],[[1]]], [[[1],[2,10],[8,10]],[[4]]]], [5,6,7,12,13,14,3,9],   [[[2],[8]], [[2,10],[1]], [[8,10],[1]], [[2,10,11],[4]], [[8,10,11],[4]]]],
                       [[[[7],[1],[0]],[[7],[1],[4]],[[7],[1],[5]],[[7],[2],[0]],[[7],[2],[4]],[[7],[2],[5]]], [[[[1],[2]],[[8]]], [[[8],[1,3],[2,3]],[[4]]]],   [5,6,7,12,13,14,9,10],  [[[1],[2]], [[1,3],[8]], [[2,3],[8]], [[1,3,11],[4]], [[2,3,11],[4]]]],
                       
                       [[[[5],[4],[0]],[[5],[4],[1]],[[5],[4],[2]],[[5],[7],[0]],[[5],[7],[1]],[[5],[7],[2]]], [[[[1],[4]],[[2]]], [[[2],[1,5],[4,5]],[[8]]]],   [9,10,11,12,13,14,3,6], [[[1],[4]], [[1,5],[2]], [[4,5],[2]], [[1,5,7],[8]], [[4,5,7],[8]]]],
                       [[[[7],[4],[0]],[[7],[4],[1]],[[7],[4],[2]],[[7],[5],[0]],[[7],[5],[1]],[[7],[5],[2]]], [[[[1],[2]],[[4]]], [[[4],[1,3],[2,3]],[[8]]]],   [9,10,11,12,13,14,5,6], [[[1],[2]], [[1,3],[4]], [[2,3],[4]], [[1,3,7],[8]], [[2,3,7],[8]]]],
                       [[[[4],[5],[0]],[[4],[5],[1]],[[4],[5],[2]],[[4],[7],[0]],[[4],[7],[1]],[[4],[7],[2]]], [[[[2],[4]],[[1]]], [[[1],[2,6],[4,6]],[[8]]]],   [9,10,11,12,13,14,3,5], [[[2],[4]], [[2,6],[1]], [[4,6],[1]], [[2,6,7],[8]], [[4,6,7],[8]]]],
                       
                       [[[[7],[1]],[[7],[2]],[[7],[4]],[[7],[5]],[[0],[1]],[[0],[2]],[[0],[4]],[[0],[5]]], [[[[1],[2]],[[4,12],[8,12]]], [[[4],[8]],[[1,3],[2,3]]]], [7,11,13,14,5,6,9,10],  [[[1],[2]], [[4],[8]], [[1,3],[4,12]], [[1,3],[8,12]], [[2,3],[4,12]], [[2,3],[8,12]]]],
                       [[[[5],[0]],[[5],[2]],[[5],[4]],[[5],[7]],[[1],[0]],[[1],[2]],[[1],[4]],[[1],[7]]], [[[[1],[4]],[[2,10],[8,10]]], [[[2],[8]],[[1,5],[4,5]]]], [7,11,13,14,3,6,9,12],  [[[1],[4]], [[2],[8]], [[1,5],[2,10]], [[1,5],[8,10]], [[4,5],[2,10]], [[4,5],[8,10]]]],
                       [[[[2],[0]],[[2],[1]],[[2],[5]],[[2],[7]],[[4],[0]],[[4],[1]],[[4],[5]],[[4],[7]]], [[[[1],[8]],[[2,6],[4,6]]], [[[2],[4]],[[1,9],[8,9]]]],   [7,11,13,14,3,5,10,12], [[[1],[8]], [[2],[4]], [[1,9],[2,6]], [[1,9],[4,6]], [[8,9],[2,6]], [[8,9],[4,6]]]]]
    
    chunks_list = []
    chunks_list_motifs = []
    for topo_in in topology_files_in:
        chunks_topo = []
        chunks_topo_motifs = []
        with open(topo_in, 'rb') as f_topo_in:
            for line in f_topo_in:
                split_line = line[0:-1].split('\t')
                chunk_motif = np.array(split_line[ind_column + 1:ind_column + 1 + 16], dtype = float)
                if np.sum(chunk_motif) < min_unmasked or split_line[0] in chromosomes_to_exclude:
                    pass
                else:
                    chunks_topo.append(['_'.join(split_line[0:ind_column + 1]), np.array(split_line[ind_column + 1:ind_column + 1 + 16], dtype = float), np.zeros(16, dtype = float), np.zeros(10, dtype = float), np.zeros(15, dtype = bool)])
                    chunks_topo_motifs.append(chunk_motif)
        chunks_topo_motifs = np.array(chunks_topo_motifs)
        chunks_list.append(chunks_topo)
        chunks_list_motifs.append(chunks_topo_motifs)
    #Fill out the mismatches
    for D_idx in range(len(chunks_list)):
        num_chunks = len(chunks_list[D_idx])
        for chunk_idx in range(num_chunks):
            bases = float(np.sum(chunks_list[D_idx][chunk_idx][1]))
            chunks_list[D_idx][chunk_idx][2] = chunks_list[D_idx][chunk_idx][1] / bases
    #Fill out the pairwise differences
    for D_idx in range(len(chunks_list)):
        num_chunks = len(chunks_list[D_idx])
        for chunk_idx in range(num_chunks):
            for pair_idx in range(len(pairwise_keys)):
                pair = pairwise_keys[pair_idx]
                tot_mismatch = 0.0
                if pair[1] == 'anc':
                    for mis_idx in range(16):
                        if mismatch_keys[mis_idx][pair[0]] == '1':
                            tot_mismatch += chunks_list[D_idx][chunk_idx][2][mis_idx]
                    chunks_list[D_idx][chunk_idx][3][pair_idx] = tot_mismatch
                else:
                    for mis_idx in range(16):
                        if mismatch_keys[mis_idx][pair[0]] != mismatch_keys[mis_idx][pair[1]]:
                            tot_mismatch += chunks_list[D_idx][chunk_idx][2][mis_idx]
                    chunks_list[D_idx][chunk_idx][3][pair_idx] = tot_mismatch
    topo_conditions_to_test_list = [['mismatch_order'], ['branch_length_order'], ['inconsistency_threshold'], ['branch_length_equality']]
    
    for topo_conditions_to_test in topo_conditions_to_test_list:
        #Fill out the topology consistencies!
        print topo_conditions_to_test
        for D_idx in range(len(chunks_list)):
            print '\t'.join(['set'] + topo_order + ["", "total"])
            num_chunks = len(chunks_list[D_idx])
            for chunk_idx in range(num_chunks):
                for topo_idx in range(len(topo_conditions)):
                    topo_condition = topo_conditions[topo_idx]
                    keep_topo = True
                    if 'mismatch_order' in topo_conditions_to_test:
                        for condition in topo_condition[0]:
                            vals = []
                            for case in condition:
                                vals.append(np.average([chunks_list[D_idx][chunk_idx][3][i] for i in case]))
                            if np.sum(np.diff(vals) > 0) != len(vals) - 1:
                                keep_topo = False
                                break
                    if 'branch_length_order' in topo_conditions_to_test and keep_topo == True:
                        #Continue checking to confirm the branch lengths are consistent
                        for condition in topo_condition[1]:
                            vals = []
                            for case in condition:
                                #Average these vals
                                to_av = []
                                for to_sum in case:
                                    to_av.append(np.sum([chunks_list[D_idx][chunk_idx][2][i] for i in to_sum]))
                                vals.append(np.average(to_av))
                            if np.sum(np.diff(vals) > 0) != len(vals) - 1:
                                keep_topo = False
                                break
                    if 'inconsistency_threshold' in topo_conditions_to_test and keep_topo == True:
                        inconsistent_motifs = np.array(topo_condition[2])
                        inconsistent = np.sum(chunks_list[D_idx][chunk_idx][2][inconsistent_motifs])
                        consistent = np.sum(chunks_list[D_idx][chunk_idx][2][np.setdiff1d(np.arange(1,15), inconsistent_motifs)])
                        if inconsistent > consistent * inconsistency_threshold:
                            keep_topo = False
                    if 'branch_length_equality' in topo_conditions_to_test and keep_topo == True:
                        for condition in topo_condition[3]:
                            vals = []
                            for to_sum in condition:
                                vals.append(np.sum([chunks_list[D_idx][chunk_idx][2][i] for i in to_sum]))
                            if abs(np.log(vals[0] / float(vals[1]))) > abs(np.log(branch_length_equality_factor)):
                                keep_topo = False
                                break
                            
                    if keep_topo == True:
                        #Topology is valid!
                        chunks_list[D_idx][chunk_idx][4][topo_idx] = True
                    else:
                        chunks_list[D_idx][chunk_idx][4][topo_idx] = False
        print "Topology counts, allowing overlapping cases (%s):" %(','.join(topo_conditions_to_test))
        topo_consistency_list = []
        for i in range(len(chunks_list)):
            topo_consistency = np.array([chunks_list[i][j][4] for j in range(len(chunks_list[i]))])
            topo_consistency_list.append(topo_consistency)
            print '\t'.join([set_names[i]] + ['%.3f' %(j) for j in np.sum(topo_consistency,0) / float(np.sum(np.sum(topo_consistency,1)>0))]) + '\t||\t%d' %(np.sum(np.sum(topo_consistency,1)>0))
        print "Topology counts, not allowing overlapping cases (%s):" %(','.join(topo_conditions_to_test))
        for i in range(len(chunks_list)):
            print '\t'.join([set_names[i]] + ['%.3f' %(j) for j in np.sum(topo_consistency_list[i][np.sum(topo_consistency_list[i],1)==1],0) / float(len(topo_consistency_list[i][np.sum(topo_consistency_list[i],1)==1]))]) + '\t||\t%d' %(len(topo_consistency_list[i][np.sum(topo_consistency_list[i],1)==1]))
        print '\t'.join(["set"] + ["0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111","","total"])
        
        for tops in topology_focus:
            print "Median x1e4 divergence of with evidence of topology %s" %(','.join([topo_order[top] for top in tops]))
            for i in range(len(chunks_list)):
                valid_mask = np.sum([topo_consistency_list[i][::,top] == True for top in tops], axis = 0, dtype = bool)
                motifs = np.array([chunks_list[i][j][2] for j in range(len(chunks_list[i]))])[valid_mask]
                print '\t'.join([set_names[i]] + ['%.2f' %(j * 1e4) for j in np.median(motifs, 0)]) + '\t||\t%.2f|%d' %(np.sum(valid_mask) / float(len(valid_mask)), np.sum(valid_mask))
        #print '\t'.join(['set'] + ["ChuDen","ChuNea","ChuHum","ChuAnc","DenNea","DenHum","DenAnc","NeaHum","NeaAnc","HumAnc","","total"])
        print '\n'
        
    return chunks_list, topo_consistency_list
