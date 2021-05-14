import sys
import time

from tqdm import tqdm
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
        nodes, elements = import_2dm_mesh(path_mesh, path_dat)
    elif path_mesh.endswith(".dat"):
        nodes, elements = import_uro_mesh(path_mesh, path_dat)
    else:
        sys.exit("Mesh format not supported.")
    print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    # print("- Reading dat...")
    # t0 = time.time()
    # print("  {}".format(path_dat))
    # import_2dm_dat(path_dat, nodes)
    # print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    mesh = Mesh(nodes, elements)

    # print("- SKIPPING steps, which are necessary for dissolving.")
    do_dissolve = True
    if do_dissolve:
        print("- Finding each nodes connected elements...")
        t0 = time.time()
        find_nodes_connected_elements(mesh)
        print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

        print("- Finding each elements neighbours...")
        t0 = time.time()
        set_element_neighbour_ids(mesh)
        print("  - Finished after {:.2f} seconds.".format(time.time() - t0))

    print("Import finished after {:.2f} seconds.".format(time.time() -
                                                         t0_import))
    return mesh


def import_2dm_mesh(path_mesh, path_dat):
    nodes, elements, elements_quad = [], [], []
    for line in open(path_mesh, encoding="latin1"):
        if line.startswith("ND "):
            nodes.append([float(x) for x in line.split()[2:]])
        elif line.startswith("E3T "):
            elements.append([int(nid)-1 for nid in line.split()[2:5]])
        elif line.startswith("E4Q "):
            nids = [int(nid)-1 for nid in line.split()[2:6]]
            elements_quad.append(nids)
    nodes = np.array(nodes)

    import_2dm_dat(path_dat, nodes)

    if len(elements_quad) > 0:
        split_method = "split_into_four"
        if split_method == "split_with_min_dz":
            print("- Splitting quads (connecting nodes with smallest dz)...")
            for nids in elements_quad:
                qnodes = nodes[nids, :]
                if abs(qnodes[0, 2] - qnodes[2, 2]) < abs(qnodes[1, 2] - qnodes[3, 2]):
                    elements.append([nids[0], nids[1], nids[2]])
                    elements.append([nids[2], nids[3], nids[0]])
                else:
                    elements.append([nids[0], nids[1], nids[3]])
                    elements.append([nids[1], nids[2], nids[3]])
        elif split_method == "split_into_four":
            print("- Splitting each quad into four triangles...")
            new_nid = nodes.shape[0]
            new_nodes = [] 
            for nids in elements_quad:
                qnodes = nodes[nids, :]
                new_nodes.append(np.mean(qnodes, axis=0))
                elements.append([nids[0], nids[1], int(new_nid)])
                elements.append([nids[1], nids[2], int(new_nid)])
                elements.append([nids[2], nids[3], int(new_nid)])
                elements.append([nids[3], nids[0], int(new_nid)])
                new_nid += 1
            nodes = np.vstack((nodes, np.array(new_nodes)))

    elements = np.array(elements)

    return nodes, elements


def import_2dm_dat(path_dat, nodes):
    # TODO: Check which dat format it is
    lines = open(path_dat, encoding="latin1").readlines()

    data = []
    for line in lines:
        try:
            data.append(float(line))
        except:
            print(f"Skipping line: {line.strip()}")
    data = np.array(data)
    nodes[:, 2] = data


def import_2dm_dat_old(path_dat, nodes):
    # TODO: Check which dat format it is
    lines = open(path_dat, encoding="latin1").readlines()

    n_nodes = nodes.shape[0]
    data = np.array([float(s) for s in lines[4:(n_nodes+4)]])
    nodes[:, 2] = data


def import_uro_mesh(path_mesh, path_dat):
    # print("- Lese UnRunOff-Mesh...")
    lines = [line for line in open(path_mesh).readlines()
             if not line.startswith('C')]

    n_bndnodes = int(lines[0])
    n_nodes = int(lines[1])
    n_elmts = int(lines[2 + n_nodes + n_bndnodes])

    i_nodes = [2, 2 + n_nodes + n_bndnodes]
    i_elmts = [2 + n_nodes + n_bndnodes + 1, 2 + n_nodes + n_bndnodes + 1 +
               n_elmts]

    def line_to_node(line):
        return [float(x) for x in line.split()[1:]]

    def line_to_elmt(line):
        return [int(nid) for nid in line.split()[:3]]

    nodes = [line_to_node(line) for line in lines[i_nodes[0]:i_nodes[1]]]
    nodes_array = np.array(nodes)

    elements = [line_to_elmt(line) for line in lines[i_elmts[0]:i_elmts[1]]]
    import_2dm_dat(path_dat, nodes_array)

    return nodes_array, np.array(elements)


def find_nodes_connected_elements(mesh):
    mesh.nodes_conn_element_ids = [[] for x in range(mesh.nodes.shape[0])]

    for eid, nids in enumerate(tqdm(mesh.elements)):
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

    for eid, nids in enumerate(tqdm(mesh.elements)):
        temp_element_neigh_ids = []
        for nid in nids:
            temp_element_neigh_ids += mesh.nodes_conn_element_ids[nid]

        for neid in temp_element_neigh_ids:
            if eid == neid:
                continue
            if len(set(nids) & set(mesh.elements[neid])) == 2:
                mesh.element_neigh_ids[eid].append(neid)
    for eid in range(len(mesh.element_neigh_ids)):
        mesh.element_neigh_ids[eid] = list(set(mesh.element_neigh_ids[eid]))



def write_element_contour_polygons(contour_polygons_dict, path_cpolys):
    print("\nWriting each elements contour polygon to shapefile...")
    print(path_cpolys)
    w = shapefile.Writer(path_cpolys)
    w.field("ID", 'N')
    w.field("Contour", 'C')

    for c_pair, contour_polygons in tqdm(contour_polygons_dict.items()):
        for i, cpoly in enumerate(contour_polygons):
            if len(cpoly) == 0:
                continue
            elif len(cpoly) > 2:
                pnts = [(x, y, z) for x, y, z in cpoly]
                pnts.append(pnts[0])
                w.polyz([pnts])
                w.record(i, c_pair)
            else:
                print("Contour Polygon {} of element {} ".format(c_pair, i) +
                      "only has {} points".format(len(cpoly)))

    w.close()
