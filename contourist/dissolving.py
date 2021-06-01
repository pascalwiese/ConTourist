import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from numpy.lib.function_base import insert
import shapefile

sys.setrecursionlimit(10000)


def dissolve_all_contour_polygons(element_contour_polygons_dict, mesh, params, do_plot=False):
    print("\nDissolving all contour polygons...")

    debug_points = []
    debug_diss_points = []
    dissolve_polygons_dict = {}
    for c_lo, c_hi in zip(params['contours'][:-1], params['contours'][1:]):
        c_key = "{}-{}".format(c_lo, c_hi)
        print("Starting with contour {}".format(c_key))
        dissolve_polygons_dict[c_key] = []
        element_contour_polygons = element_contour_polygons_dict[c_key]
        is_elmt_dissolved = [False for i in range(len(mesh.elements))]

        for eid in range(len(is_elmt_dissolved)):
            if len(element_contour_polygons[eid]) == 0:
                is_elmt_dissolved[eid] = True

        print("- Calculating the dissolve order...")
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

            dissolve_order = [eid]
            dissolve_order = get_dissolve_order_recursive(dissolve_order, is_elmt_dissolved, eid, element_contour_polygons, mesh)

            all_dissolve_orders.append(dissolve_order)
        print("  -> {} areas found.".format(len(all_dissolve_orders)))

        print("- Dissolving all areas to polygons...")
        for dissolve_order in all_dissolve_orders:
            polygon = dissolve_polygons(element_contour_polygons, dissolve_order, c_key, do_plot=do_plot)
        

        with open("./testdata/diss_order_{}.txt".format(c_key), "w") as fout:
            fout.write("id ord eid x y z\n")
            for i, dissolve_order in enumerate(all_dissolve_orders):
                for ord, eid in enumerate(dissolve_order):
                    x, y, z = element_to_point(mesh.elements[eid], mesh.nodes)
                    fout.write(f"{i} {ord} {eid} {x} {y} {z}\n")
                      

def dissolve_polygons(polygons, dissolve_order, c_key, do_plot=False):
    polygon = polygons[dissolve_order[0]]

    for eid in dissolve_order[1:]:
        # print(polygons[eid])
        polygon = dissolve(polygon, polygons[eid], do_plot=do_plot)

    return polygon 


# Dissolving two polygons that share at least two vertices
def dissolve(poly1, poly2, do_plot=False):
    # if is_polygon_ccw(poly1) != is_polygon_ccw(poly2):
    #     poly2 = poly2[::-1]
    
    print("poly1:", poly1)
    print("poly2:", poly2)
    poly1_in_poly2 = [None for i in range(len(poly1))]
    poly2_in_poly1 = [None for i in range(len(poly2))]

    for i2, pnt2 in enumerate(poly2):
        # print(pnt2)
        for i1, pnt1 in enumerate(poly1):
            # print(pnt1)
            if pnt2[0] == pnt1[0] and pnt2[1] == pnt1[1]:
                poly1_in_poly2[i1] = i2
                poly2_in_poly1[i2] = i1
    print("poly1_in_poly2:", poly1_in_poly2)
    print("poly2_in_poly1:", poly2_in_poly1)
    insert_idx = 0
    for i1, i2 in zip(poly1_in_poly2[:-1], poly1_in_poly2[1:]):
        insert_idx += 1
        if i1 is not None and i2 is not None:
            break

    print("insert_idx:", insert_idx)
    print()

    if True:
        fig, ax = plt.subplots()
        patches = []
        np_poly = np.array([[pnt[0], pnt[1]] for pnt in poly1])
        ax.scatter(np_poly[:, 0], np_poly[:, 1])
        patches.append(Polygon(np_poly, False))

        np_poly = np.array([[pnt[0], pnt[1]] for pnt in poly2])
        ax.scatter(np_poly[:, 0], np_poly[:, 1])
        patches.append(Polygon(np_poly, False))

        p = PatchCollection(patches, alpha=0.5)
        colors = 100 * np.random.rand(len(patches))
        p.set_array(colors)
        ax.add_collection(p)
        fig.colorbar(p, ax=ax)
        ax.axis("equal")
        plt.show()

    # dissolve
    if None in poly2_in_poly1:
        poly = poly1
        poly2 = [pnt for pnt, inside in zip(poly2, poly2_in_poly1) if inside is None]
        poly = poly[:insert_idx] + poly2 + poly[insert_idx:]
    else:  # all points already appear in poly
        pt1 = poly1[:min(poly2_in_poly1)+1]
        pt2 = poly1[max(poly2_in_poly1):]
        # print("pt1:", pt1)
        # print("pt2:", pt2)
        poly = pt1 + pt2



    return poly


# This is only robust for convex polygons. The contour polygons for each element are convex though
def is_polygon_ccw(poly):
    a, b, c = poly[:3]
    det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
    if det > 0:
        return True
    else:
        return False
    

def get_dissolve_order_recursive(dissolve_order, is_elmt_dissolved, eid, element_contour_polygons, mesh):
    is_elmt_dissolved[eid] = True
    neids = [neid for neid in mesh.element_neigh_ids[eid] 
             if len(element_contour_polygons[neid]) > 0 and 
             not is_elmt_dissolved[neid]]
    dissolve_order += neids

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
    