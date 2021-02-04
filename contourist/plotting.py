# import matplotlib.pyplot as plt


# def plot_single_element_contour(nodes, c_poly):
#     fig, ax = plt.subplots()

#     nodes_xs = [nd[0] for nd in nodes] + [nodes[0][0]]
#     nodes_ys = [nd[1] for nd in nodes] + [nodes[0][1]]
#     ax.plot(nodes_xs, nodes_ys, color='black')
#     for x, y, z in nodes:
#         ax.text(x, y, z)

#     poly_xs = [pnt[0] for pnt in c_poly] + [c_poly[0][0]]
#     poly_ys = [pnt[1] for pnt in c_poly] + [c_poly[0][1]]
#     ax.plot(poly_xs, poly_ys, color='red')
#     for x, y, z in c_poly:
#         ax.text(x, y, z)

#     plt.show()
