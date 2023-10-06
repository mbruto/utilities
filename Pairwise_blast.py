#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

import os
import argparse
import sys
import re
import numpy
import scipy.cluster.hierarchy as sciclus
import matplotlib.pyplot as plt
import pylab
import subprocess

#-------------------

def prepare_data(sub_csv_file) :

    hash_spec = {}
    
    with open(sub_csv_file) as file :
        for line in file :
            if line[0] == '#':
                continue
            
            temp_list = line.split(",")
            species_name = temp_list[0].split(" ")[0]+"_"+temp_list[0].split(" ")[1]
            species_name = re.sub("\"", "", species_name)
            #print(species_name)
            if species_name in hash_spec :
                hash_spec[species_name].append(temp_list[5][1:-1])
            else :
                hash_spec[species_name] = [temp_list[5][1:-1]]
    
    return hash_spec
    
#-------------------

def filter_data(sub_in_dict, sub_in_fasta) :

    id_list = []
    to_del_list = []

    with open(sub_in_fasta) as file :
        for line in file :
            if line[0] == '>' :
                id_list.append(line[1:].split("_")[0]+"_"+line[1:].split("_")[1])
    
    for sp in sub_in_dict.keys() :
        gca_to_remove = []
        for gca in sub_in_dict[sp] :
            if gca not in id_list :
               gca_to_remove.append(gca)
        for gca in gca_to_remove :
            sub_in_dict[sp].remove(gca)
        
        if len(sub_in_dict[sp]) == 0 :
            to_del_list.append(sp)
            
    for d in to_del_list :
        sub_in_dict.pop(d, None)
    
    return sub_in_dict

#-------------------

def get_blast_res(sub_blast_file) :

    sub_blast_dict = {}
    
    with open(sub_blast_file) as b_file :
        for line in b_file :
            if not line.rstrip() :
                continue
            
            temp_list = line.split("\t")
            
            sub_gca = temp_list[0].split("_")[0]+"_"+temp_list[0].split("_")[1]
            
            if sub_gca in sub_blast_dict :
                if temp_list[0] in sub_blast_dict[sub_gca] :
                    sub_blast_dict[sub_gca][temp_list[0]][temp_list[1]] = int(temp_list[4]) + int(temp_list[5])
                else :
                    sub_blast_dict[sub_gca][temp_list[0]] = { temp_list[1] : int(temp_list[4]) + int(temp_list[5]) }
            else :
                sub_blast_dict[sub_gca] = { temp_list[0] : {temp_list[1] : int(temp_list[4]) + int(temp_list[5]) }}
    
    return sub_blast_dict


#-------------------

def ANI_matrix(sub_in_fasta, sub_in_dict) :


    gen_list = list(sub_in_dict.keys())
    
    #Counting the number of genomes
    genome_num = len((sub_in_dict.keys()))
    
    #Create the matrix containing the ANI results
    ani_array = numpy.zeros((genome_num, genome_num))

    #Launch fastANI and index the results
    x_index = -1
    
    cmd_blast = "blastn -task blastn -query "+sub_in_fasta+" -subject "+sub_in_fasta+" -max_target_seqs 1500 -qcov_hsp_perc 80 -outfmt 6 -out temp_result.csv"
    os.system(cmd_blast)
    
    blast_dict = get_blast_res("temp_result.csv")
    #print(blast_dict['GCA_001267155.1']['GCA_001267155.1_0']['GCA_002237575.1_0'])
    #exit()
    
    for i in range(0, len(gen_list)) :
       
        species_1 = gen_list[i]
        print(species_1)
    
        for j in range(0, len(gen_list)) :
       
            species_2 = gen_list[j]
            count = 0
            total_polym = 0
            
            #print("\t"+species_2)
        
            for strain_i in sub_in_dict[species_1] :
                
                for strain_j in sub_in_dict[species_2] :
                    
                    if strain_i == strain_j :
                        total_polym += 0
                        count += 1
                        continue
                    
                    for rrna_i in blast_dict[strain_i].keys() :
                        
                        for rrna_j in blast_dict[strain_i][rrna_i] :
                            #print(rrna_j+"\t"+strain_j)
                            if rrna_j.startswith(strain_j) :
                                count += 1
                                total_polym += blast_dict[strain_i][rrna_i][rrna_j]

                if count == 0 :
                    ani_array[i, j] = 500.0
                    #ani_array[j, i] = 0.0
                else :
                    ani_array[i, j] = total_polym / count
                    #ani_array[j, i] = total_polym / count
                
    
    print(ani_array)

    return ani_array, gen_list
#-------------------

def create_figs(sub_ani_mat, sub_name_list) :

    clust_name_list = []
    
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(12,12))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6], frame_on=False)
    Y = sciclus.linkage(sub_ani_mat, method='complete')
    Z1 = sciclus.dendrogram(Y, orientation='left', distance_sort="ascending")
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    #Reorder the lit of names according to the clustering
    for i in Z1['leaves'] :
       clust_name_list.append(sub_name_list[i])
    
    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2], frame_on=False)
    Y = sciclus.linkage(sub_ani_mat, method='complete')
    Z2 = sciclus.dendrogram(Y, distance_sort = "descending")
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    sub_ani_mat = sub_ani_mat[idx1,:]
    sub_ani_mat = sub_ani_mat[:,idx2]
    im = axmatrix.matshow(sub_ani_mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    
    #Set names
    axmatrix.set_xticks(range(len(clust_name_list)))
    axmatrix.set_xticklabels(reversed(clust_name_list), minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()

    pylab.xticks(rotation=-90, fontsize=5)

    axmatrix.set_yticks(range(len(clust_name_list)))
    axmatrix.set_yticklabels(clust_name_list, minor=False, fontsize=5)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    
    # Plot colorbar. left, bottom, width, height
    axcolor = fig.add_axes([0.09,0.71,0.1,0.2])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('dendrogram.png')
    
    return clust_name_list

#-------------------

def output_matrix(sub_ani_mat, sub_name_list, sub_order_list) :

    OUT = open("ANI_matrix.csv", "w")
    sub_order_list.reverse()
    
    for n in sub_order_list :
        OUT.write("\t"+str(n))
    OUT.write("\n")
    
    for i in range(0,len(sub_order_list)) :
        OUT.write(sub_order_list[i])
        index_i = sub_name_list.index(sub_order_list[i])
        for j in range(0,len(sub_name_list)) :
            index_j = sub_name_list.index(sub_order_list[j])
            OUT.write("\t"+str(sub_ani_mat[index_i][index_j]))
        OUT.write("\n")
    
    OUT.close()

#-------------------

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Create a matrix of pairwise sequence identity')
    parser.add_argument('-i', '--input', required=True, help="Input fasta file")
    #parser.add_argument('-c', '--csv', required=True, help="CSV NCBI file")
    parser.add_argument('-o', '--out', required=False, help="Output files name")
    args = parser.parse_args()
    
    
    taxo_dict = prepare_data(args.csv)
    taxo_dict = filter_data(taxo_dict, args.input)
    ANI_mat, name_list = ANI_matrix(args.input, taxo_dict)
    order_name_list = create_figs(ANI_mat, name_list)
    output_matrix(ANI_mat, name_list, order_name_list)
