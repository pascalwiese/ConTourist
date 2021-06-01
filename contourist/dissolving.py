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
            polygon = dissolve_polygons(element_contour_polygons, dissolve_order, c_key)
        

        with open("./testdata/diss_order_{}.txt".format(c_key), "w") as fout:
            fout.write("id ord eid x y z\n")
            for i, dissolve_order in enumerate(all_dissolve_orders):
                for ord, eid in enumerate(dissolve_order):
                    x, y, z = element_to_point(mesh.elements[eid], mesh.nodes)
                    fout.write(f"{i} {ord} {eid} {x} {y} {z}\n")
                      

def dissolve_polygons(polygons, dissolve_order, c_key):
    polygon = polygons[dissolve_order[0]]

    for eid in dissolve_order[1:]:
        # print(polygons[eid])
        polygon = dissolve(polygon, polygons[eid])

    return polygon 


# Dissolving two polygons that share at least two vertices
def dissolve(poly1, poly2):
    # if is_polygon_ccw(poly1) != is_polygon_ccw(poly2):
    #     poly2 = poly2[::-1]
    
    print("poly1:", poly1)
    print("poly2:", poly2)
    poly1_in_poly2 = [0 for i in range(len(poly1))]
    poly2_in_poly1 = [0 for i in range(len(poly2))]

    for i2, pnt2 in enumerate(poly2):
        # print(pnt2)
        for i1, pnt1 in enumerate(poly1):
            # print(pnt1)
            if pnt2[0] == pnt1[0] and pnt2[1] == pnt1[1]:
                poly1_in_poly2[i1] = 1
                poly2_in_poly1[i2] = 1
    print("poly1:", poly1_in_poly2)
    print("poly2:", poly2_in_poly1)
    print()

    poly = [poly1]
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
    