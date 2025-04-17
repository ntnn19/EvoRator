import sys
from .paths import path_to_cdhit,path_to_mafft
import numpy as np
import os
from ..preprocessing.sequence_utils import load_FASTA, num2seq
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from datetime import datetime
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import itertools
def write_alignment(all_sequences, location):
    with open(location, 'w') as f:
        for k, sequence in enumerate(all_sequences):
            f.write('>%s\n' % k)
            f.write(sequence + '\n')
    return location



def apply_cdhit(all_sequences, threshold=0.95, length_cutoff=0.8, name=None, cdhit=path_to_cdhit):
    if name is None:
        timestamp = str(datetime.now()).replace(':', '_').replace(' ', '_')
        name = 'sequences_%s' % timestamp

    cdhit_input = 'MSA_%s.fasta' % name
    cdhit_output = 'representatives_%s.fasta' % name
    cdhit_output_clusters = cdhit_output + '.clstr'
    write_alignment(all_sequences, cdhit_input)
    if threshold > 0.7:
        n = 5
    elif threshold > 0.6:
        n = 4
    elif threshold > 0.5:
        n = 3
    elif threshold > 0.4:
        n = 2
    else:
        print('Threshold %s too low, aborting' % threshold)
        return

    command = '%s -i %s -o %s -c %s -n %s -T 0 -M 0 -s %s -g 1 -G 0 -aL %s -aS %s -l 5' % (
        cdhit, cdhit_input, cdhit_output, threshold, n, length_cutoff,length_cutoff,length_cutoff)
    os.system(command)

    B = len(all_sequences)
    labels = np.zeros(B, dtype=np.int)
    nClusters = 0
    nSequences_per_cluster = []
    all_representatives = []

    cluster = -1
    nSequences = 0
    with open(cdhit_output_clusters, 'r') as f:
        for line in f:
            if line[0] == '>':  # New cluster.
                if cluster >= 0:  # Add result of previous cluster.
                    nSequences_per_cluster.append(nSequences)
                    nClusters += 1
                cluster = int(line[9:-1])
                nSequences = 0
            elif len(line) > 1:
                sequence = int(line.split('>')[1].split('...')[0])
                labels[sequence] = cluster
                nSequences += 1
                if '*' in line:
                    all_representatives.append(sequence)

        nSequences_per_cluster.append(nSequences)
        nClusters += 1

    all_representatives = np.array(all_representatives)
    nSequences_per_cluster = np.array(nSequences_per_cluster)

    os.system('rm %s' % cdhit_input)
    os.system('rm %s' % cdhit_output)
    os.system('rm %s' % cdhit_output_clusters)
    return labels, nClusters, nSequences_per_cluster, all_representatives


def apply_cdhit_safe(all_sequences, threshold=0.95,threshold_safe=None, length_cutoff=0.8, name=None, cdhit=path_to_cdhit):
    if threshold_safe is None:
        if threshold >=0.8:
            threshold_safe = 0.7
        else:
            threshold_safe = threshold - 0.1
    labels, nClusters, nSequences_per_cluster, all_representatives = apply_cdhit(all_sequences, threshold=threshold,
                                                                                 length_cutoff=length_cutoff, name=name,
                                                                                 cdhit=cdhit)
    all_representative_sequences = np.array([all_sequences[x] for x in all_representatives])
    labels_bis, nClusters_bis, nSequences_per_cluster_bis, all_representatives_bis = apply_cdhit(
        all_representative_sequences, threshold=threshold_safe, length_cutoff=length_cutoff, name=name, cdhit=cdhit)

    total_comparisons = 0
    clusters_to_merge = []

    for cluster in range(nClusters_bis):
        cluster_size = nSequences_per_cluster_bis[cluster]
        original_clusters = np.nonzero(labels_bis == cluster)[0]
        if cluster_size > 1:
            comparisons = (cluster_size * (cluster_size - 1)) // 2
            print('Performing %s pairwise comparison' % comparisons)
            total_comparisons += comparisons
            labels_, nClusters_, clusterSize_, _ = pairwise_cluster_sequences(
                all_representative_sequences[original_clusters], threshold=threshold, length_cutoff=length_cutoff)
            if nClusters_ < cluster_size:
                print('Found %s missing mergers' % (cluster_size - nClusters_))
                for n in range(nClusters_):
                    if clusterSize_[n] > 1:
                        clusters_to_merge.append(original_clusters[(labels_ == n)])

    labels_final = labels.copy()

    for cluster_to_merge in clusters_to_merge:
        merger = cluster_to_merge[0]
        for mergee in cluster_to_merge[1:]:
            labels_final[labels_final == mergee] = merger
    unique_values, labels_final, nSequences_per_cluster_final = np.unique(labels_final, return_inverse=True,
                                                                          return_counts=True)
    nClusters_final = len(unique_values)
    all_representatives_final = []
    for n in range(nClusters_final):
        subset = np.nonzero(labels_final == n)[0]
        representative = subset[np.argmax([len(x) for x in all_sequences[subset]])]
        all_representatives_final.append(representative)
    all_representatives_final = np.array(all_representatives_final)

    print('First CD-HIT found %s clusters, performed %s pairwise comparisons and refined to %s clusters' % (
    nClusters, total_comparisons, nClusters_final))
    return labels_final, nClusters_final, nSequences_per_cluster_final, all_representatives_final

def apply_mafft(sequences, mafft=path_to_mafft,go_penalty=1.53,
    ge_penalty= 0.0,name=None, numeric=False, return_index=False,high_accuracy=False):
    if name is None:
        name = '%.10f' % np.random.rand()
    input_file = 'tmp_%s_unaligned.fasta' % name
    output_file = 'tmp_%s_aligned.fasta' % name
    instruction_file = 'tmp_%s.sh' % name
    write_alignment(sequences, input_file)
    if high_accuracy:
        command = '%s  --amino --localpair --maxiterate 1000 --op %s --ep %s %s > %s' % (mafft, go_penalty,ge_penalty,input_file, output_file)
    else:
        command = '%s  --amino --auto --op %s --ep %s %s > %s' % (mafft, go_penalty,ge_penalty,input_file, output_file)
    print(command)
    with open(instruction_file, 'w') as f:
        f.write(command)
    os.system('sh %s' % instruction_file)
    alignment = load_FASTA(
        output_file, drop_duplicates=False)[0]
    if return_index:
        is_gap = alignment == 20
        index = np.cumsum(1 - is_gap, axis=1) - 1
        index[is_gap] = -1

    if not numeric:
        alignment = num2seq(alignment)
    os.system('rm %s' % input_file)
    os.system('rm %s' % output_file)
    os.system('rm %s' % instruction_file)

    if return_index:
        return alignment, index
    else:
        return alignment


def align_and_map(seq1, seq2, return_similarity=False,local=False):
    # a = pairwise2.align.globalxx(seq1, seq2)[0]
#     a = pairwise2.align.globalxs(seq1,seq2,-4,-0.1)[0]
    if local:
        a = pairwise2.align.localds(seq1, seq2, blosum62, -10, -0.5, one_alignment_only=True)[0]
    else:
        a = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5,one_alignment_only=True)[0]
    lalign = a[-1]
    sites1 = []
    sites2 = []
    count1 = -1
    count2 = -1
    for l in range(lalign):
        aa1 = a[0][l]
        aa2 = a[1][l]
        if aa1 != '-':
            count1 += 1
        if aa2 != '-':
            count2 += 1
        if aa1 == aa2:
            sites1.append(count1)
            sites2.append(count2)
    sites1 = np.array(sites1)
    sites2 = np.array(sites2)
    if return_similarity:
        similarity = len(sites1)/ min([len(seq1), len(seq2)])
        return sites1, sites2, similarity
#         return sites1, sites2, a[2]
    return sites1, sites2

def pairwise_cluster_sequences(all_sequences,threshold=0.95,length_cutoff=0.8):
    nSequences = len(all_sequences)
    pairwise_seqIds = np.eye(nSequences)
    pairwise_coverage = np.eye(nSequences)
    lengths = np.array([len(sequence) for sequence in all_sequences])
    for n1 in range(nSequences):
        for n2 in range(n1+1,nSequences):
            sites1, sites2, similarity = align_and_map(all_sequences[n1],all_sequences[n2],return_similarity=True)
            pairwise_seqIds[n1,n2] = similarity
            pairwise_coverage[n1,n2] = len(sites1)/len(all_sequences[n1])
            pairwise_coverage[n2, n1] = len(sites2) / len(all_sequences[n2])

    pairwise_seqIds += pairwise_seqIds.T
    graph = (pairwise_seqIds>= threshold) & (pairwise_coverage>=length_cutoff) & (pairwise_coverage.T>=length_cutoff)
    nClusters, labels = connected_components(csgraph=csr_matrix(graph), return_labels=True)
    nSequences_per_cluster = np.zeros(nClusters)
    all_representatives = []

    for n in range(nClusters):
        indices = np.nonzero(labels == n)[0]
        nSequences_per_cluster[n] = len(indices)
        representative = indices[np.argmax(lengths[indices])]
        all_representatives.append(all_sequences[representative])
    return labels,nClusters,nSequences_per_cluster,all_representatives


def apply_cdhit_double_safe(all_sequences, all_organisms, all_caths, threshold=0.95, threshold_safe=0.7,
                            length_cutoff=0.8, name=None, cdhit=path_to_cdhit,
                            watch_keywords=[]):
    def construct_list_to_compare(all_organisms, all_caths,
                                  watch_keywords=[]):
        def simplify_organism(x):
            x = x.split(' (')[0]
            x = [y for y in x.split(' ') if not '/' in y]
            x = ' '.join(x[:2])
            x = x.lower()
            return x

        nSequences = len(all_organisms)

        all_organisms = np.array(all_organisms)
        organism_simplified = np.array(list(map(simplify_organism, all_organisms)))

        organisms = np.unique(organism_simplified)
        if len(watch_keywords) > 0:
            watched_organisms = [organism for organism in organisms if any([x in organism for x in watch_keywords])]
        else:
            watched_organisms = organisms

        filter_watched_organisms = np.array([organism in watched_organisms for organism in organism_simplified])

        same_organism = np.array([
            [organism_simplified[i] == organism_simplified[j] for j in range(nSequences)] for i in range(nSequences)
        ])
        same_organism = same_organism & filter_watched_organisms[np.newaxis] & filter_watched_organisms[:, np.newaxis]
        all_caths = np.array(all_caths)
        same_cath = np.array(
            [[(all_caths[i] in all_caths[j]) | (all_caths[j] in all_caths[i]) for j in range(nSequences)] for i in
             range(nSequences)])
        same_cath = same_cath & (all_caths != '')[np.newaxis] & (all_caths != '')[:, np.newaxis]

        compare = same_cath | same_organism
        Is, Js = np.nonzero(compare)
        list_to_compare = np.array([(i, j) for i, j in zip(Is, Js) if j > i])
        return list_to_compare

    labels, nClusters, nSequences_per_cluster, all_representatives = apply_cdhit_safe(all_sequences,
                                                                                      threshold=threshold,
                                                                                      threshold_safe=threshold_safe,
                                                                                      length_cutoff=length_cutoff,
                                                                                      name=name, cdhit=cdhit)
    lengths = np.array([len(x) for x in all_sequences])
    list_to_compare = construct_list_to_compare(all_organisms, all_caths, watch_keywords=watch_keywords)
    print('%s additional comparisons to be performed' % len(list_to_compare) )
    cluster_graph = np.eye(nClusters, dtype=np.bool)

    count_pairwise = 0
    count_cluster = dict([('%s_%s'%(i,j) , 0) for i,j in itertools.product(np.arange(nClusters),np.arange(nClusters) )])

    for k,pair in enumerate(list_to_compare):
        i,j = pair
        if (labels[i] != labels[j]) & (cluster_graph[labels[i], labels[j]] == 0) & (count_cluster['%s_%s'%(labels[i],labels[j])] < 20):
            count_pairwise += 1
            count_cluster['%s_%s' % (labels[i], labels[j])] +=1
            if count_pairwise % 100 ==1:
                print('Performed %s additional pairwise comparisons (%.2f)' % (count_pairwise, k/len(list_to_compare) ) )

            sites1, sites2, similarity = align_and_map(all_sequences[i], all_sequences[j], return_similarity=True)
            coverage1 = len(sites1) / len(all_sequences[i])
            coverage2 = len(sites2) / len(all_sequences[j])
            if (similarity >= threshold) & (coverage1 >= length_cutoff) & (coverage2 >= length_cutoff):
                cluster_graph[labels[i], labels[j]] = 1
                cluster_graph[labels[j], labels[i]] = 1
                print('Merge found!')
    # cluster_graph += cluster_graph.T

    nClusters_final, labels_ = connected_components(csgraph=csr_matrix(cluster_graph), return_labels=True)
    labels_final = labels_[labels]
    nSequences_per_cluster_final = np.zeros(nClusters_final)

    all_representatives_final = []

    for n in range(nClusters_final):
        indices = np.nonzero(labels_final == n)[0]
        nSequences_per_cluster_final[n] = len(indices)
        representative = indices[np.argmax(lengths[indices])]
        all_representatives_final.append(representative)
    all_representatives_final = np.array(all_representatives_final)

    print('Second CD-HIT found %s clusters, performed %s pairwise comparisons and refined to %s clusters' % (
    nClusters, count_pairwise, nClusters_final))
    return labels_final, nClusters_final, nSequences_per_cluster_final, all_representatives_final
