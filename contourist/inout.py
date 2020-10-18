import sys
import time

import numpy as np
import shapefile

from .classes import Mesh


def import_mesh_data(path_mesh, path_dat):
    print("\nImporting mesh and data...")
    t0_import = time.time()

    print("- Reading mesh...")
    t0 = time.time()
    print("  {}".format(path_mesh))
    # TODO: Check if mesh is renumbered
    if path_mesh.endswith(".2dm"):
        nodes, elements = import_2dm_mesh(path_mesh)
    else:
        sys.exit("Mesh format not supported.")
    print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    print("- Reading dat...")
    t0 = time.time()
    print("  {}".format(path_dat))
    import_2dm_dat(path_dat, nodes)
    print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    mesh = Mesh(nodes, elements)

    print("- Finding each nodes connected elements...")
    t0 = time.time()
    find_nodes_connected_elements(mesh)
    print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    print("- Finding each elements neighbours...")
    t0 = time.time()
    set_element_neighbour_ids(mesh)
    print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    print("- Import finished after {:.2f} seconds.".format(time.time() - 
                                                           t0_import))
    return mesh


def import_2dm_mesh(path_mesh):
    nodes, elements = [], []
    for line in open(path_mesh, encoding="latin1"):
        if line.startswith("ND "):
            nodes.append([float(x) for x in line.split()[2:]])
        elif line.startswith("E3T "):
            elements.append([int(nid)-1 for nid in line.split()[2:5]])
        elif line.startswith("E4Q "):
            nids = [int(nid)-1 for nid in line.split()[2:6]]
            elements.append([nids[0], nids[1], nids[2]])
            elements.append([nids[2], nids[3], nids[0]])

    nodes = np.array(nodes)
    elements = np.array(elements)

    return nodes, elements


def import_2dm_dat(path_dat, nodes):
    # TODO: Check which dat format it is
    lines = open(path_dat, encoding="latin1").readlines()

    n_nodes = nodes.shape[0]
    data = np.array([float(s) for s in lines[4:(n_nodes+4)]])
    nodes[:, 2] = data


def find_nodes_connected_elements(mesh):
    mesh.nodes_conn_element_ids = [[] for x in range(mesh.nodes.shape[0])]

    for eid, nids in enumerate(mesh.elements):
        mesh.nodes_conn_element_ids[nids[0]].append(eid)
        mesh.nodes_conn_element_ids[nids[1]].append(eid)
        mesh.nodes_conn_element_ids[nids[2]].append(eid)

    # TODO: Check, if the following is necessary.
    #       Totally forgot, if it is... >_>
    mesh.nodes_conn_element_ids = [list(set(eids)) for eids in
                                   mesh.nodes_conn_element_ids]


# For each element the neighbouring elements are found
# Neighbouring means that they are connected by an edge
# and not only one node
def set_element_neighbour_ids(mesh):
    mesh.element_neigh_ids = [[] for x in range(mesh.elements.shape[0])]

    for eid, nids in enumerate(mesh.elements):
        temp_element_neigh_ids = []
        for nid in nids:
            temp_element_neigh_ids += mesh.nodes_conn_element_ids[nid]

        for neid in temp_element_neigh_ids:
            if eid == neid:
                continue
            if len(set(nids) & set(mesh.elements[neid])) == 2:
                mesh.element_neigh_ids[eid].append(neid)


def write_element_contour_polygons(contour_polygons, path_cpolys):
    print("\nWriting each elements contour polygon to shapefile...")
    print(path_cpolys)
    w = shapefile.Writer(path_cpolys)
    w.field("ID", 'N')

    for i, cpoly in enumerate(contour_polygons):
        w.poly([pnt for pnt in cpoly])
        w.record(i)
    w.close()
