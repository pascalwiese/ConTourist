import os
import argparse
from sys import path

import numpy as np

from contourist.inout import import_mesh_data
from contourist.inout import write_element_contour_polygons
from contourist.contouring import create_contour_polygons


def main():
    # Parsing command line arguments
    path_dict, params = parse_args()

    # Reading mesh and data
    mesh = import_mesh_data(path_dict['mesh'], path_dict['dat'])

    # Creating contours for each element
    element_contour_polygon_dict = create_contour_polygons(mesh, params,
                                                           debug=False)

    # Writing undissolved contour polygons to shapefile
    write_element_contour_polygons(element_contour_polygon_dict,
                                   path_dict['out'])


def parse_args():
    description = ("Beschreibung folgt")
    # Construct the argument parser and parse the arguments
    ap = argparse.ArgumentParser(description=description,
                                 epilog="$ConTourist",
                                 formatter_class=argparse.RawTextHelpFormatter)

    ap.add_argument("--mesh", metavar='netz.2dm', type=str,
                    required=True,
                    help="Name oder Pfad des Meshs (Hydro_AS, UnRunOff)")
    ap.add_argument("--dat", metavar='DEPTH.dat', type=str,
                    required=False, default=None,
                    help="Name oder Pfad der skalaren Ergebnisdatei (.dat)")
    ap.add_argument("--contours", metavar='0,0.5,1,999', type=str,
                    required=True, default='0,0.5,1,999',
                    help="Kontur-Klassen (mit Komma getrennt)")
    ap.add_argument("--out", metavar="contours.shp", type=str, required=False,
                    default=None, help="Optionaler Ausgabedateiname")
    ap.add_argument("--diss", metavar='false', type=str,
                    required=False, default='false',
                    help="Dissolve die Einzelkonturflaechen (true/false)")
    ap.add_argument("--plot", metavar='false', type=str,
                    required=False, default='false',
                    help="Kontroll-Plots erstellen.")
    args = vars(ap.parse_args())

    folder = '.'
    if os.path.isabs(args['mesh']):
        path_mesh = args['mesh']
    else:
        path_mesh = os.path.join(folder, args['mesh'])

    if os.path.isabs(args['dat']):
        path_dat = args['dat']
    else:
        path_dat = os.path.join(folder, args['dat'])

    if args['out']:
        path_out = os.path.join(os.path.dirname(path_dat), args['out'])
    else:
        path_out = path_dat[:-4] + "_contours.shp"

    print("\nAusgabepfad:\n{}".format(path_out))

    contours = [float(x) for x in args['contours'].split(',')]
    do_plot = True if args['plot'].upper() in ["TRUE", "JA", "J", "Y"] \
        else False
    do_diss = True if args['diss'].upper() in ["TRUE", "JA", "J", "Y"] \
        else False

    paths = {'mesh': path_mesh, 'dat': path_dat, 'out': path_out}
    params = {'contours': contours, 'plot': do_plot, 'dissolve': do_diss}

    return paths, params


if __name__ == '__main__':
    main()
