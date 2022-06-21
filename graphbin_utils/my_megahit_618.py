from Bio import SeqIO
from bidirectionalmap.bidirectionalmap import BidirectionalMap
import re
from igraph import *
if __name__ == '__main__':
        # Get original contig IDs
    #-------------------------------

    original_contigs = {}
    #contigs文件
    for index, record in enumerate(SeqIO.parse('/media/ubuntu/conda/vicent/HiCzin/my_data_617_megahit/k141.final.contigs.fa', "fasta")):
        original_contigs[record.id] = str(record.seq)

    contig_descriptions = {}

    for index, record in enumerate(SeqIO.parse('/media/ubuntu/conda/vicent/HiCzin/my_data_617_megahit/k141.final.contigs.fa', "fasta")):
        contig_descriptions[record.id] = record.description

    # Construct the assembly graph
    #-------------------------------

    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

        # Get links from .gfa file
    with open('/media/ubuntu/conda/vicent/HiCzin/my_data_617_megahit/k141.gfa') as file:

        line = file.readline()

        while line != "":

            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                start_1 = 'NODE_'
                end_1 = '_length'

                link1 = int(re.search('%s(.*)%s' % (start_1, end_1), strings[1]).group(1))

                start_2 = 'NODE_'
                end_2 = '_length'

                link2 = int(re.search('%s(.*)%s' % (start_2, end_2), strings[3]).group(1))

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                strings = line.split()

                start = 'NODE_'
                end = '_length'

                contig_num = int(re.search('%s(.*)%s' % (start, end), strings[1]).group(1))

                my_map[node_count] = int(contig_num)

                graph_contigs[contig_num] = strings[2]

                node_count += 1

            line = file.readline()

            contigs_map = my_map
        contigs_map_rev = my_map.inverse

    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Create list of edges
    edge_list = []

    for i in range(node_count):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= str(contigs_map[i])

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))
            edge = str(contigs_map_rev[link[0]])+'  ' + \
            str(contigs_map_rev[link[1]])+'   1'+"\n"
            with open("intermediate.txt", "a") as f:
                f.write(edge)


    import shutil
    readDir = "intermediate.txt"
    writeDir = "assembly_Graph_vicent618.txt"
    lines_seen = set()
    outfile = open(writeDir, "w")
    f = open(readDir, "r")
    for line in f:
        # 删除相同行
        if line not in lines_seen:
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()
    # print "success"

    # Add edges to the graph
    # print(edge_list)

    # assembly_graph.add_edges(edge_list)
    # # print(assembly_graph)
    # # assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)
    # # print(type(assembly_graph))
    # print(str(assembly_graph))

    # with open("assembly_graph_vicent.txt", "w") as f:
    #     f.write(str(assembly_graph))

