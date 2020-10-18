import numpy as np
from tqdm import tqdm

from .geometry import point_on_line_by_z
from .plotting import plot_single_element_contour


def create_contour_polygons(mesh, params, debug=False):
    print("\nCreating contours polygons for each element...")
    element_contour_polygons = {}

    for c_lo, c_hi in zip(params['contours'][:-1],
                          params['contours'][1:]):
        print("- Starting with contour {} <= contour < {}...".format(c_lo,
                                                                     c_hi))
        contour_polygons = [[] for x in range(mesh.elements.shape[0])]
       
        for eid, nids in enumerate(tqdm(mesh.elements)):
            nodes = mesh.nodes[nids, :]

            # All z values are lower than the contour bounds
            if np.sum(nodes[:, 2] < c_lo) == 3:
                continue
            # All z values are higher than the contour bounds
            elif np.sum(nodes[:, 2] > c_hi) == 3:
                continue
            # All z values are within the contour bounds
            elif np.sum(nodes[:, 2] >= c_lo) == 3 and \
                    np.sum(nodes[:, 2] <= c_hi) == 3:
                contour_polygons[eid] = nodes
            # Element must have a contour polygon
            else:
                for nd1, nd2 in zip(nodes, [nodes[1], nodes[2], nodes[0]]):
                    # Check, if node pair is outside contour bounds
                    if nd1[2] > c_hi and nd2[2] > c_hi:
                        continue
                    if nd1[2] < c_lo and nd2[2] < c_lo:
                        continue
                    
                    # Checking which nodes' z value is higher
                    nd_lo = nd1
                    nd_hi = nd2
                    is_reverse = False
                    if nd_lo[2] > nd_hi[2]:
                        nd_lo, nd_hi = nd_hi, nd_lo
                        is_reverse = True

                    # Both nodes inside contour bounds.
                    # => Trivial case
                    if nd_lo[2] >= c_lo and nd_hi[2] <= c_hi:
                        if is_reverse:
                            contour_polygons[eid] += [nd2, nd1]
                        else:
                            contour_polygons[eid] += [nd1, nd2]
                    # Otherwise:
                    # => One or two intersections possible
                    else:
                        intersections = []
                        # Intersection with lower contour value
                        if nd_lo[2] < c_lo:
                            intersections.append(point_on_line_by_z(nd_lo,
                                                                    nd_hi,
                                                                    c_lo))
                        else:
                            intersections.append(nd_lo)

                        # Intersection with higher contour value
                        if nd_hi[2] > c_hi:
                            intersections.append(point_on_line_by_z(nd_lo,
                                                                    nd_hi,
                                                                    c_hi))
                        else:
                            intersections.append(nd_hi)

                        if is_reverse:
                            contour_polygons[eid] += intersections[::-1]
                        else:
                            contour_polygons[eid] += intersections

                if debug:
                    print("\nElement index: {}".format(eid))
                    for pnt in contour_polygons[eid]:
                        print(pnt)
                    plot_single_element_contour(nodes,
                                                contour_polygons[eid])

        # Adding all contour polygons to the dictionary
        element_contour_polygons["{}-{}".format(c_lo, c_hi)] = contour_polygons

    return element_contour_polygons
