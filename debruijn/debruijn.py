#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

from random import randint
import statistics
import textwrap
import argparse
import os
import sys
import random
import networkx as nx
#import matplotlib
from operator import itemgetter
random.seed(9001)
#import matplotlib.pyplot as plt
#matplotlib.use("Agg")

__author__ = "Emeline PAPILLON"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Emeline PAPILLON"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Emeline PAPILLON"
__email__ = "emelinepap@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as files_fastq:
        for i in files_fastq:
            yield next(files_fastq).strip("\n")
            next(files_fastq)
            next(files_fastq)


def cut_kmer(read, kmer_size):
    for i in range(0, len(read)-kmer_size+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    dic = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in dic:
                dic[kmer] += 1
            else:
                dic[kmer] = 1
    return dic


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        graph.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    if (statistics.stdev(weight_avg_list)) > 0:
        max_weight_index = weight_avg_list.index(max(weight_avg_list))
        del path_list[max_weight_index]
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    elif (statistics.stdev(path_length)) > 0:
        max_length_index = path_length.index(max(path_length))
        del path_list[max_length_index]
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    else:
        rand = randint(0, len(path_list))
        del path_list[rand]
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    weight_avg_list = []
    path_length = []
    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        weight_avg_list.append(path_average_weight(graph, path))
        path_length.append(len(path))
    return select_best_path(graph, path_list, path_length, weight_avg_list)


def simplify_bubbles(graph):
    bubble = False
    for node in graph.nodes():
        list_predecessor = list(graph.predecessors(node))
        if len(list_predecessor) > 1:
            for i in list_predecessor:
                for j in list_predecessor[list_predecessor.index(i) + 1:len(list_predecessor)]:
                    ancestor_node = nx.lowest_common_ancestor(graph, i, j)
                    if ancestor_node is not None:
                        bubble = True
                        break
        if bubble:
            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    start = False
    path_list = []
    path_length = []
    weight_avg_list = []
    for node in graph.nodes():
        list_predecessor = list(graph.predecessors(node))
        if len(list_predecessor) > 1:
            for starting_nod in starting_nodes:
                if len(list(nx.all_simple_paths(graph, starting_nod, node))) != 0:
                    start = True
                    for path in nx.all_simple_paths(graph, starting_nod, node):
                        path_list.append(path)
                        weight_avg_list.append(path_average_weight(graph, path))
                        path_length.append(len(path))
    if start:
        graph = select_best_path(graph, path_list, path_length, weight_avg_list,
                                 delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph, ending_nodes):
    end = False
    path_list = []
    path_length = []
    weight_avg_list = []
    for node in graph.nodes():
        list_successors = list(graph.successors(node))
        if len(list_successors) > 1:
            for ending_nod in ending_nodes:
                if len(list(nx.all_simple_paths(graph, node, ending_nod))) != 0:
                    end = True
                    for path in nx.all_simple_paths(graph, node, ending_nod):
                        path_list.append(path)
                        weight_avg_list.append(path_average_weight(graph, path))
                        path_length.append(len(path))
    if end:
        graph = select_best_path(graph, path_list, path_length, weight_avg_list,
                                 delete_entry_node=False, delete_sink_node=True)
    return graph

def get_starting_nodes(graph):
    entry_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            entry_nodes.append(node)
    return entry_nodes

def get_sink_nodes(graph):
    sink_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    list_contigs = []
    for n_entry in starting_nodes:
        for n_sink in ending_nodes:
            if nx.has_path(graph, n_entry, n_sink):
                for path in nx.all_simple_paths(graph, n_entry, n_sink):
                    contig = path[0]
                    for node in path[1:]:
                        contig += node[-1]
                    list_contigs.append([contig, len(contig)])
    return list_contigs

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as file:
        for i in range(len(contigs_list)):
            file.write(">contig_{} len={}\n{}\n".format(i, contigs_list[i][1],
                                                        textwrap.fill(contigs_list[i][0],
                                                                      width=80)))


def draw_graph(graph, graphimg_file):
    """Draw the graph"""
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    """
    #Test for a little file
    for read in read_fastq(args.fastq_file):
        #print(read)
        cut_kmer(read, args.kmer_size)
    """

    #Lecture du fichier et construction du graphe
    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))

    #Résolution des bulles
    simplify_bubbles(graph)

    #Résolution des portes d'entrée et de sorties
    solve_entry_tips(graph, get_starting_nodes(graph))
    solve_out_tips(graph, get_sink_nodes(graph))

    #Ecriture du/des contigs
    save_contigs(get_contigs(graph, get_starting_nodes(graph),
                             get_sink_nodes(graph)), args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
