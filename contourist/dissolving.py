import sys

import numpy as np
import matplotlib.pyplot as plt
import shapefile

sys.setrecursionlimit(10000)


def dissolve_all_contour_polygons(element_contour_polygons_dict, mesh, params, do_plot=False):
    print("\nDissolving all contour polygons...")

    debug_points = []
    debug_diss_points = []
    dissolve_polygons_dict = {}
    for c_lo, c_hi in zip(params['contours'][:-1], params['contours'][1:]):
        c_key = "{}-{}".format(c_lo, c_hi)
        print("- Starting with contour {}".format(c_key))
        dissolve_polygons_dict[c_key] = []
        element_contour_polygons = element_contour_polygons_dict[c_key]
        is_elmt_dissolved = [False for i in range(len(mesh.elements))]

        for eid in range(len(is_elmt_dissolved)):
            if len(element_contour_polygons[eid]) == 0:
                is_elmt_dissolved[eid] = True

        all_dissolve_orders = []
        for eid in range(len(is_elmt_dissolved)):
            # Skip if element was dissolved
            if is_elmt_dissolved[eid]:
                continue

            # Skip if element has no contour polygon 
            if len(element_contour_polygons[eid]) == 0:
                is_elmt_dissolved[eid] = True
                continue
            
            is_elmt_dissolved[eid] = True

            print("\n", eid)
            neids = [neid for neid in mesh.element_neigh_ids[eid] 
                     if len(element_contour_polygons[neid]) > 0 and 
                     not is_elmt_dissolved[neid]]

            dissolve_order = [eid]
            dissolve_order = get_dissolve_order_recursive(dissolve_order, is_elmt_dissolved, eid, element_contour_polygons, mesh)

            print(dissolve_order)
            all_dissolve_orders.append(dissolve_order)

        with open("./testdata/diss_order_{}.txt".format(c_key), "w") as fout:
            fout.write("id ord eid x y z\n")
            for i, dissolve_order in enumerate(all_dissolve_orders):
                for ord, eid in enumerate(dissolve_order):
                    x, y, z = element_to_point(mesh.elements[eid], mesh.nodes)
                    fout.write(f"{i} {ord} {eid} {x} {y} {z}\n")
                      


def get_dissolve_order_recursive(dissolve_order, is_elmt_dissolved, eid, element_contour_polygons, mesh):
    is_elmt_dissolved[eid] = True
    neids = [neid for neid in mesh.element_neigh_ids[eid] 
             if len(element_contour_polygons[neid]) > 0 and 
             not is_elmt_dissolved[neid]]
    
    dissolve_order += neids
    
    # print("- {} - {} - {}".format(eid, neids, dissolve_order))
    # a = input()
    for neid in neids:
        is_elmt_dissolved[neid] = True
        
    for neid in neids:    
        get_dissolve_order_recursive(dissolve_order, is_elmt_dissolved, neid, element_contour_polygons, mesh)

    return dissolve_order
    

def dissolve_with_neighbours(polygon, eid, neids, element_contour_polygons, c_key, mesh, do_plot=False):
    print("n_neids: {}".format(len(neids)))
    if len(neids) == 0:
        print("-- finished")
        return polygon


def recursive_dissolving(polygon, eid, element_contour_polygons, is_elmt_dissolved, c_key, mesh, do_plot=False):
    is_elmt_dissolved[eid] = True
    polygon.append(element_to_point(mesh.elements[eid], mesh.nodes))

    dissolve_neighbour_eids = [neid for neid in mesh.element_neigh_ids[eid] 
                               if len(element_contour_polygons[neid]) > 0 and 
                               not is_elmt_dissolved[neid]]

    # print("\nElement: ", eid)
    # print("Neighbours: ", mesh.element_neigh_ids[eid])
    # print("Dissolve Neighbours: ", dissolve_neighbour_eids)
    
    if False:
        for pnt in polygon:
            plt.scatter(pnt[0], pnt[1], color="yellow", s=120)

        xs, ys = element_to_lines(mesh.elements[eid], mesh.nodes)
        plt.plot(xs, ys, color="black")
        centroid = element_to_point(mesh.elements[eid], mesh.nodes)
        plt.scatter(centroid[0], centroid[1], color="red")
        # plt.text(str(eid+1), centroid[0], centroid[1])

        for neid in mesh.element_neigh_ids[eid]:
            centroid = element_to_point(mesh.elements[neid], mesh.nodes)
            plt.scatter(centroid[0], centroid[1], color="black", s=10)
            xs, ys = element_to_lines(mesh.elements[neid], mesh.nodes)
            plt.plot(xs, ys, color="black")
            # plt.text(str(eid+1), centroid[0], centroid[1])

        for neid in dissolve_neighbour_eids:
            centroid = element_to_point(mesh.elements[neid], mesh.nodes)
            plt.scatter(centroid[0], centroid[1], color="blue")

        plt.axis("equal")
        plt.show()


    # If no more neighbours can be dissolved
    if len(dissolve_neighbour_eids) == 0:
        return polygon

    # First do dissolving
    for neid in dissolve_neighbour_eids:
        # diss_polygon.append(element_to_point(mesh.elements[neid], mesh.nodes))
        # polygon.append(element_to_point(mesh.elements[neid], mesh.nodes))
        is_elmt_dissolved[neid] = True

    # Then dissolve the neighbours neighbours
    for neid in dissolve_neighbour_eids:
        polygon = recursive_dissolving(polygon, neid, element_contour_polygons, is_elmt_dissolved, c_key, mesh, do_plot=True)
    
    return polygon


def element_to_point(element, nodes):
    elmt_nodes = [nodes[nid] for nid in element]
    x = sum([n[0] for n in elmt_nodes]) / len(elmt_nodes)
    y = sum([n[1] for n in elmt_nodes]) / len(elmt_nodes)
    z = sum([n[2] for n in elmt_nodes]) / len(elmt_nodes)
    return [x, y, z]


def element_to_lines(element, nodes):
    xs = [nodes[nid][0] for nid in element]
    ys = [nodes[nid][1] for nid in element]
    xs.append(xs[0])
    ys.append(ys[0])
    return (xs, ys)


def write_debug_points_to_shape(all_points, name="debug"):
    w = shapefile.Writer("./testdata/{}_lines.shp".format(name))
    w.field('LID', 'N')

    for i, points in enumerate(all_points):
        w.line([points])
        w.record(i)
    w.close()

    w2 = shapefile.Writer("./testdata/{}_points.shp".format(name))
    w2.field('LID')
    w2.field('PID')

    for i, points in enumerate(all_points):
        for pi, pnt in enumerate(points):
            w2.point(pnt[0], pnt[1])
            w2.record(i, pi)
    w2.close()
    