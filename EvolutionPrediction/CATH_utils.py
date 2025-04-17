import numpy as np

# https://www.cathdb.info/wiki?id=data:index
# https://filezilla-project.org
import os
import sys
laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True
if laptop:
    location_cath = "cath/cath_b.20201021.all"
    location_cath_labels = "cath/cath_b.names.20201021"
else:
    location_cath = "/bioseq/evorator/EvolutionPrediction/cath/cath_b.20201021.all"
    location_cath_labels = "/bioseq/evorator/EvolutionPrediction/cath/cath_b.names.20201021"



def load_cath_database(location_cath=location_cath):
    result_dict = {}
    with open(location_cath, 'r') as f:
        for line in f:
            line = line[:-1]  # Remove \n.
            parsed_line = line.split(' ')
            pdb = parsed_line[0][:4]
            chain = parsed_line[0][4]
            identifier = '%s_%s' % (pdb, chain)
            version = parsed_line[1]
            cath = parsed_line[2]
            location_line = ' '.join(parsed_line[3:]).split(', ')
            locations = []
            for word in location_line:
                word = word.split(':')[0]
                end = word.split('-')[-1]
                start = '-'.join(word.split('-')[:-1])
                locations.append((start, end))
            all_cath = [(cath, location) for location in locations]
            if identifier in result_dict.keys():
                result_dict[identifier] += all_cath
            else:
                result_dict[identifier] = all_cath
    for key, item in result_dict.items():  # Reorder.
        if len(item) > 1:
            start = []
            for dom in item:
                try:
                    start.append(int(dom[1][0]))
                except:
                    start.append(int(dom[1][0][:-1]))
            order = np.argsort(start)
            item_sorted = [item[order_] for order_ in order]
            result_dict[key] = item_sorted
    return result_dict


def load_cath_labels(location_cath_labels):
    labels_dict = {}
    with open(location_cath_labels, 'r',encoding='utf8') as f:
        for line in f:
            line = line[:-1]
            parsed_line = line.split(' ')
            cath_id = parsed_line[0]
            cath_label = ' '.join(parsed_line[1:])
            labels_dict[cath_id] = cath_label
    return labels_dict


cath_labels = load_cath_labels(location_cath_labels)
cath_database = load_cath_database(location_cath)

def get_cath_class(cath_ids, cath_labels=cath_labels):
    islist = isinstance(cath_ids, list)
    if not islist:
        cath_ids = [cath_ids]
    all_labels = []
    for cath_id in cath_ids:
        label = []
        for level in range(4):
            partial_id = '.'.join(cath_id.split('.')[:level + 1])
            partial_label = cath_labels[partial_id]
            label.append(partial_label)
        all_labels.append(label)
    if not islist:
        return all_labels[0]
    else:
        return all_labels
