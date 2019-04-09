### This function counts and writes out the occurence of different mutation motifs in a
### comparison between an target individual and other comparison individuals for blocks
### specific in a BED file. The comparison individuals are often Denisovan, Neanderhtal
### and an outgroup human individual.

import time
import os
import tempfile
import subprocess
import numpy as np
import copy
import re
import gzip


def calculateMotifs_BEDonBCF(bed = None, ind_name = None, outfile = 'calculateTopology_%s.out.bed', full_bcf_file = False, comparison_list = ['DenisovaPinky', 'AltaiNea', 'LP6005441-DNA_G10'], bcf_base = '/path_to_bcf_directory/bcf_filename_chromosome%d.bcf.bgzf', random_chunks = 'disabled', bedtools_dir = './bedtools2/bin/', tmp_folder = './tmp/', inds_are_chroms = True):
    """
    This is to count up the number of times positions have each mutation motifs in a comparison of a specified set of human and archaic genomes.

    For example, we might be interested in the 16 possible ancestry-aware combinations in [chunk, deni, nean, human] order.
    Human should be an individual (often an African but potentially a European or other individual, depending on the use).
    Based on the number of times we observe each mutation motif we can approximately assess support for different coalescent topologies for the blocks of interest.

    bed is the set of blocks to determine topology for
    ind_name is the target individual name in the BCF
    comparison_list is the list of individuals in the BCF file that are used to define a mutation motif. E.g. ['DENI', 'NEAN', 'HUMAN_OUTGROUP'].
    
    If inds_are_chroms is True then this assumes the BED file was calculated on the same phased data as the bcf, and that archaic haplotypes have been predicted in a phase aware way.
    In this case, the ind_name should be given as IND_1 for the first (phased) chromosome copy and IND_2 for the second, where IND is the name of the target individual in the BCF.

    If inds_are_chroms is False then this assumes the BED file makes unphased predictions of archaic haplotypes.
    Thus, it isn't meaningful to differentiate chromosome copies. Mutation motif counts will be far harder to interpret as topologies in this circumstance as the two chromosome copies will often have different coalescent topologies.
    Now ind_name can be given as IND, where IND is the name of the target individual in the BCF.
    Two files are written, the motif counts on the first (_1) and second (_2) chromosomes, based on the same input BED file.

    full_bcf_file says whether or not I'm calculating topologies on a BCF file that keeps invariable sites. In this case, I don't filter down to variable SNPs when trimming out individuals not relevant to the calculation. This can be useful e.g. if some sites are called and some are not in the BCF file and we want to keep track of how many sites were called, but the BCF files will be larger and the calculation slower.
    Note that there is an extra column when full_bcf_file is True, The final column now counts how many bases were included in the motif counting for the block. Note that e.g. alignability masking especially will often mean that not all SNPs in a block are called.
    
    I save the chr, chunk start and end position, and the number of mutation motifs of each type observed.
    It is assumed that I will use the AA here; I can fold it down later.

    Make repeated calls to singleAssess_mismatchTopology to retrieve the topology of each chunk, by default based on a chunk, Deni, Nean, San comparison.
    Alternatively, if mismatches are already computed, computed_file = FILENAME simply read this from a corresponding file FILENAME %(ind)
    If using S* (inds_are_chroms = False) then the assessment will be made for each chromosome and I have the option of discarding the more-human chromosome later.
    In each case, save out a list of chrom, start, end, mismatches.
    The key order is 0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
    
    """
    start_time = time.time()
    print "Beginning mismatch toplogy assessment of ind %s at time" %(ind_name), start_time
    BEDTOOLS_DIR = bedtools_dir
    TMP_FOLDER = tmp_folder
    if not os.path.exists(tmp_folder):
        os.mkdir(tmp_folder)
    ind_offset = [0,1] if inds_are_chroms == False else [0] if ind_name[-2:] == '_1' else [1] if ind_name[-2:] == '_2' else np.nan #offset for the vcf processing
    #Convert ind into a temp file
    command_convert = "cat %s" %(bed)
    if random_chunks == 'shuffle':
        #In this case, shuffle the actual intersecting chunks after conversion. Bit slower but more accurate!
        command_convert = command_convert + PIPE + BEDTOOLS_DIR + "bedtools shuffle -i stdin -g %s -excl %s -noOverlapping" %(GENOME_FILE, HG19_MASK)
    f_chunks = tempfile.NamedTemporaryFile(mode = 'w+b', suffix = '.bed', dir = TMP_FOLDER)
    bedtools_conv_ab = subprocess.Popen(args = command_convert, stdout = f_chunks, stderr = f_chunks, shell = True)
    bedtools_conv_ab.communicate(input=None)
    f_chunks.flush()
    f_chunks.seek(0)
    
    intersect_chunk_dict = [{i+1:[] for i in range(22)}, {i+1:[] for i in range(22)}]
    #Replace chr[xx] with [xx] in the BED file.
    for line in f_chunks:
        split_line = line[0:-1].split('\t')
        intersect_chunk_dict[0][int(split_line[0].split('chr')[-1])].append([int(split_line[1]), int(split_line[2])])
        intersect_chunk_dict[1][int(split_line[0].split('chr')[-1])].append([int(split_line[1]), int(split_line[2])])
    samples_to_keep = [ind_name[0:-2] if inds_are_chroms == True else ind_name]
    for i in comparison_list:
        if '=' not in i and i != 'REF':
            samples_to_keep.append(i)
    samples_to_keep = ','.join(samples_to_keep)

    #Prep dict
    num_topology = [[0], [1]]
    for comp in range(len(comparison_list)):
        new_topology = []
        for top in num_topology:
            new_topology.append(top + [0])
            new_topology.append(top + [1])
        num_topology = copy.copy(new_topology)
    num_topology = {''.join(['%d' %(i) for i in top]):[0.0,0.0] for top in num_topology}
    num_topology_keys = np.sort(num_topology.keys())
    length_mismatch = [[],[]] #_0 and _1; only one is used if not S*
    for chrom in range(1,23):
        print "Cutting down chromosome %d BCF to archaic chunks, target inds are %s and %s" %(chrom, ind_name, samples_to_keep)
        #intersect_chunk_dict[chrom] = np.array(intersect_chunk_dict[chrom], dtype = float)
        #Use bcftools to create a temporary VCF file, containing variant bi-allelic snps only
        f_intersectvcf = tempfile.NamedTemporaryFile(mode = 'w+b', dir = TMP_FOLDER)
        #Command converted to bcftools. Note that more recent bcftools might allow trimming by genotype = save Python steps.
        if full_bcf_file == False:
            command_intersectbcf = 'bcftools view %s --regions-file %s --samples %s --min-ac=1:minor -Ov' %(bcf_base %(chrom), f_chunks.name, samples_to_keep)
        else:
            command_intersectbcf = 'bcftools view %s --regions-file %s --samples %s -Ov' %(bcf_base %(chrom), f_chunks.name, samples_to_keep)
        vcftools_process_maxvcf = subprocess.Popen(args = command_intersectbcf, stdout = f_intersectvcf, stderr = subprocess.PIPE, shell = True)
        vcftools_process_maxvcf.communicate(input=None)
        f_intersectvcf.flush()
        f_intersectvcf.seek(0)
        
        #For each SNP, assess which category it is in for the haplotype.
        #If it is different, check if it is in the current chunk. If not, move the chunk forward until it is. Add the the total mismatches. NB that I ignore transition/transversion.
        curr_chunk = 0
        num_snps = 0
        num_topology = {top:[0.0,0.0] for top in num_topology_keys}
        num_topology['Obs'] = 0
        
        skipped = 0
        for line in f_intersectvcf:
            snp_skip = 0
            if line[0:2] == '##':
                pass
            elif line[0] == '#': #headers
                split_line = np.array(line[1:-1].split('\t'))
                target_ind = ind_name[0:-2] if inds_are_chroms == True else ind_name
                ind_idx = np.where(split_line == target_ind)[0][0]
                anc_idx = np.where(split_line == 'INFO')[0][0]
                comparative_idx_list = []
                for i in comparison_list:
                    comparative_idx_list.append(np.where(split_line == i)[0][0])
                #print split_line
                #print comparative_idx_list
                        
            else: #snp
                num_snps += 1
                split_line = line[0:-1].split('\t')
                #if split_line[0] == 'chr4' and int(split_line[1]) < 2684307:
                #    print '\t'.join(split_line)
                anc_allele = split_line[anc_idx].split('AA=')[1][0] #This is ATGC, comparative_idx is INFO field
                ref_allele = split_line[3] #VCF specification
                alt_allele = split_line[4] #VCF specification
                anc_allele = 0 if anc_allele.upper() == ref_allele.upper() else 1 if anc_allele.upper() == alt_allele.upper() else '.'
                if anc_allele == '.' or line.count('\t./.') > 0 or line.count('\t-1/-1') > 0 or '.' in split_line[ind_idx] or len(alt_allele) > 1: #one of anc, target or comparative inds is missing. NOTE that the bcf and vcf conversions appear to change missing alleles to -1s?
                    snp_skip = 1
                else:
                    #if '0' in split_line[ind_idx] and '1' in split_line[ind_idx]:
                    #    print split_line
                    hap_types = [[], []]
                    for comp_IDx in xrange(len(comparison_list)):
                        comparative_idx = comparative_idx_list[comp_IDx]
                        #print archaic_list[archaic_IDx], comparative_idx
                        comparative_alleles = re.split('[/|]', split_line[comparative_idx])
                        #NOTE that I require just use one hap here. May create issues with hets, homs and phasing. But it's probably OK.
                        if comparative_alleles[0] == '.' or comparative_alleles[1] == '.':
                            #It is possible that one of the comparisons is '.'. Can they be '-1'?
                            #Missing information - skip the line.
                            #Note that anc, ref, nean and deni will have different amounts of missing information; but this is true whatever really...
                            snp_skip = 1
                            break
                        elif comparative_alleles[0] == '-1' or comparative_alleles[1] == '-1':
                            print "BUG CATCHER!"
                        else:
                            comparative_alleles = np.array(comparative_alleles, dtype = int)
                            # I want to process both with 1/2 weighting if not phased
                            # And one with 100% weighting if phased
                            
                            if '/' not in split_line[ind_idx]:
                                #Should be phased
                                for hap in ind_offset:
                                    if len(hap_types[hap]) == 0:
                                        ind_allele = int(split_line[ind_idx].split('|')[hap])
                                        ind_allele = ind_allele if anc_allele == 0 else abs(ind_allele - 1)
                                        hap_types[hap] = [[ind_allele]]
                                    # This is assuming that the comparative allele to consider is hap 1
                                    comp_allele = comparative_alleles[0] if anc_allele == 0 else abs(comparative_alleles[0] - 1)
                                    hap_types[hap][0].append(comp_allele)
                            else:
                                # '/' is in split_line[ind_idx]
                                # Add each possible comparative allele to the mix
                                for hap in ind_offset:
                                    if len(hap_types[hap]) == 0:
                                        ind_alleles = np.array(re.split('[/|]', split_line[ind_idx]), dtype = int)
                                        ind_alleles = ind_alleles if anc_allele == 0 else np.abs(ind_alleles - 1)
                                        for ind_allele in ind_alleles:
                                            hap_types[hap].append([ind_allele])
                                    comparative_alleles = comparative_alleles if anc_allele == 0 else np.abs(comparative_alleles - 1)
                                    new_hap_types = []
                                    for possibility in hap_types[hap]:
                                        new_hap_types.append(possibility + [comparative_alleles[0]])
                                        new_hap_types.append(possibility + [comparative_alleles[1]])
                                    hap_types[hap] = copy.copy(new_hap_types)
                                    #if '0' in split_line[ind_idx] and '1' in split_line[ind_idx]:
                                    #    print hap_types
                    
                    if snp_skip == 1:
                        skipped += 1
                    else:
                        num_topology['Obs'] += 1
                        #if '0' in split_line[ind_idx] and '1' in split_line[ind_idx]:
                        #    print num_topology
                        if full_bcf_file == False:
                            for hap in ind_offset:
                                num_topology[''.join(['%d' %(i) for i in hap_types[hap][0]])][hap] += 1
                        else:
                            # Average the possibilities
                            for hap in ind_offset:
                                num_possibilities = len(hap_types[hap])
                                for possibility in hap_types[hap]:
                                    num_topology[''.join(['%d' %(i) for i in possibility])][hap] += (1.0 / num_possibilities)
                        #if split_line[0] == 'chr4' and int(split_line[1]) < 2684307:
                        #if '0' in split_line[ind_idx] and '1' in split_line[ind_idx]:
                        #    print num_topology
                        pos = int(split_line[1])
                        write_now = False
                        while pos > intersect_chunk_dict[0][chrom][curr_chunk][1]:
                            if write_now == True:
                                #We've already moved on, write that 0s were observed
                                for hap in ind_offset:
                                    if full_bcf_file == False:
                                        length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk - 1][0], intersect_chunk_dict[0][chrom][curr_chunk - 1][1]] + [0 for top in num_topology_keys])
                                    else:
                                        length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk - 1][0], intersect_chunk_dict[0][chrom][curr_chunk - 1][1]] + [0 for top in num_topology_keys] + [0])
                            curr_chunk += 1
                            write_now = True
                        #if split_line[0] == 'chr4' and int(split_line[1]) < 2704307:
                        #    print pos, intersect_chunk_dict[0][chrom][curr_chunk][1], write_now
                        if write_now == True:
                            #print "writing"
                            for hap in ind_offset:
                                if full_bcf_file == False:
                                    length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk - 1][0], intersect_chunk_dict[0][chrom][curr_chunk - 1][1]] + [num_topology[key][hap] for key in num_topology_keys])
                                else:
                                    length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk - 1][0], intersect_chunk_dict[0][chrom][curr_chunk - 1][1]] + [num_topology[key][hap] for key in num_topology_keys] + [num_topology['Obs']])
                            #print num_topology_keys
                            #print length_mismatch[0][-1]
                            num_topology = {top:[0,0] for top in num_topology_keys}
                            num_topology['Obs'] = 0
        #Last one. If there are any unwritten chunks in the chromosome, write them.
        #This is, first the one that wans't written, and then any others.
        if len(intersect_chunk_dict[0][chrom]) > 0:
            for remaining_chunk in range(len(intersect_chunk_dict[0][chrom]) - curr_chunk):
                for hap in ind_offset:
                    if full_bcf_file == False:
                        length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk + remaining_chunk][0], intersect_chunk_dict[0][chrom][curr_chunk + remaining_chunk][1]] + [num_topology[key][hap] for key in num_topology_keys])
                    else:
                        length_mismatch[hap].append([chrom, intersect_chunk_dict[0][chrom][curr_chunk + remaining_chunk][0], intersect_chunk_dict[0][chrom][curr_chunk + remaining_chunk][1]] + [num_topology[key][hap] for key in num_topology_keys] + [num_topology['Obs']])
                num_topology = {top:[0,0] for top in num_topology_keys}
                num_topology['Obs'] = 0
        f_intersectvcf.close()
    f_chunks.close()
    print "Completed assessment of ind %s in time" %(ind_name), time.time() - start_time
    keys_to_write = num_topology_keys if full_bcf_file == False else np.concatenate((num_topology_keys, np.array(['Obs'])))

    str_format = '%d' if full_bcf_file == False else '%.2f'
    num_to_write = 2**(len(comparison_list) + 1) if full_bcf_file == False else (2**(len(comparison_list) + 1)) + 1
    open_out = gzip.open if outfile[-3:] == '.gz' else open
    
    #Write out
    if inds_are_chroms == True:
        f_out = open_out(outfile %(ind_name), 'wb')
        f_out.close()
    else:
        f_out = open_out(outfile %(ind_name + '_1'), 'wb')
        f_out.close()
        f_out = open_out(outfile %(ind_name + '_2'), 'wb')
        f_out.close()
    #singleAssess_mismatchTopology(bed = assess_f, ind_name = ind, full_bcf_file = full_bcf_file, comparison_list = comparison_list, bcf_base = bcf_base, random_chunks = random_chunks, inds_are_chroms = inds_are_chroms)
    ind_result = [length_mismatch, keys_to_write]
    if inds_are_chroms == True:
        with open_out(outfile %(ind_name), 'wb') as f_out:
            assert (int(ind_name[-1]) == 1 or int(ind_name[-1]) == 2)
            for chunk in ind_result[0][int(ind_name[-1]) - 1]:
                line_to_write = ['%d' %(i) for i in chunk[0:3]] + [str_format %(chunk[3+i]) for i in range(num_to_write)]
                f_out.write(','.join(line_to_write) + '\n')
            f_out.close()
    else:
        with open_out(outfile %(ind_name + '_1'), 'wb') as f_out:
            for chunk in ind_result[0][0]:
                line_to_write = ['%d' %(i) for i in chunk[0:3]] + [str_format %(chunk[3+i]) for i in range(num_to_write)]
                f_out.write(','.join(line_to_write) + '\n')
            f_out.close()
        with open_out(outfile %(ind_name + '_2'), 'wb') as f_out:
            for chunk in ind_result[0][1]:
                line_to_write = ['%d' %(i) for i in chunk[0:3]] + [str_format %(chunk[3+i]) for i in range(num_to_write)]
                f_out.write(','.join(line_to_write) + '\n')
            f_out.close()

    return None

