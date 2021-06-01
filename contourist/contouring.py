import numpy as np
from tqdm import tqdm

from .geometry import point_on_line_by_z
# from .plotting import plot_single_element_contour


def create_contour_polygons(mesh, params, debug=False):
    print("\nCreating contours polygons for each element...")
    element_contour_polygons = {}

    for c_lo, c_hi in zip(params['contours'][:-1],
                          params['contours'][1:]):
        print("- Starting with contour {} <= contour < {}...".format(c_lo,
                                                                     c_hi))
        contour_polygons = [[] for x in range(mesh.elements.shape[0])]

        # TODO: Put in several different functions
        for eid, nids in enumerate(tqdm(mesh.elements)):
            nodes = mesh.nodes[nids, :]

            # All z values are lower than the contour bounds
            if np.sum(nodes[:, 2] <= c_lo) == 3:
                continue
            # All z values are higher than the contour bounds
            elif np.sum(nodes[:, 2] >= c_hi) == 3:
                continue
            # All z values are within the contour bounds
            elif np.sum(nodes[:, 2] >= c_lo) == 3 and \
                    np.sum(nodes[:, 2] <= c_hi) == 3:
                contour_polygons[eid] = [[n[0], n[1], n[2]] for n in nodes]
            # Element must have a contour polygon
            else:
                for nd1, nd2 in zip(nodes, [nodes[1], nodes[2], nodes[0]]):
                    # Check, if node pair is outside contour bounds
                    if nd1[2] > c_hi and nd2[2] > c_hi:
                        continue
                    if nd1[2] < c_lo and nd2[2] < c_lo:
                        continue

                    # Checking which nodes' z value is higher
                    nd_lo = [nd1[0], nd1[1], nd1[2]]
                    nd_hi = [nd2[0], nd2[1], nd2[2]]
                    is_reverse = False
                    if nd_lo[2] > nd_hi[2]:
                        nd_lo, nd_hi = nd_hi, nd_lo
                        is_reverse = True

                    # Both nodes inside contour bounds.
                    # => Trivial case. Nodes will be added with previous or
                    #    next edge.
                    if nd_lo[2] >= c_lo and nd_hi[2] <= c_hi:
                        continue

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
                            intersections = intersections[::-1]

                        contour_polygons[eid] += intersections

                if debug:
                    print("\nElement index: {}".format(eid))
                    # for pnt in contour_polygons[eid]:
                    #     print(pnt)
                    # plot_single_element_contour(nodes,
                    #                             contour_polygons[eid])

        # Deleting dublicate points
        for eid, points in enumerate(contour_polygons):
            if len(points) == 0:
                continue

            points_cleaned = [points[0]]
            for pnt in points[1:]:
                if pnt[0] != points_cleaned[-1][0] and pnt[1] != points_cleaned[-1][1]:
                    points_cleaned.append(pnt)
            if points_cleaned[0][0] == points_cleaned[-1][0] and points_cleaned[0][1] == points_cleaned[-1][1]:
                points_cleaned = points_cleaned[1:]
            contour_polygons[eid] = points_cleaned



        # Adding all contour polygons to the dictionary
        element_contour_polygons["{}-{}".format(c_lo, c_hi)] = contour_polygons
        print("  -> {} contour polygons found within {} and {}".format(len(contour_polygons), c_lo, c_hi))

    return element_contour_polygons

