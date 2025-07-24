"""
    Plot all cells on ax using patch collections for speed
"""
# Standard
from matplotlib.collections import PatchCollection, LineCollection
import numpy as np

# Custom
from RodShapedBacteria import RodShapedBacterium
from ChainingRodShapedBacteria import ChainingRodShapedBacterium

# Third party
from descartes import PolygonPatch
from shapely.geometry import Point

def add_polygon_to_plot(ax, cell, fc, ec):
    cell.addElementToPlot(ax, colour=fc, ec=ec)

def chaining_status_to_color(cell):
    """
    Returns red if the cell is not chained (no upper or lower link),
    and green if it is chained (has either upper or lower link).
    """
    upper_link = cell.upper_link
    lower_link = cell.lower_link
    radius=cell.radius
    if (upper_link is None or np.isnan(upper_link) or upper_link == "None") and (lower_link is None or np.isnan(lower_link) or lower_link == "None") and radius==0.5:
        return "#00ffff"  # Blue for Pseudomonas
    if (upper_link is None or np.isnan(upper_link) or upper_link == "None") and (lower_link is None or np.isnan(lower_link) or lower_link == "None") and radius==2.05983:
        return "#9e003a"
    if (upper_link is None or np.isnan(upper_link) or upper_link == "None") and (lower_link is None or np.isnan(lower_link) or lower_link == "None") and radius==1.00855:
        return "#9e003a"
    else:
        return (1, 0, 0, 1)  # red for Candida

    #colour the ids that cause crashing
yellow_cells_ids = [ 2331 , 2309 ] 
# def chaining_status_to_color(cell):
#     """
#     Returns yellow if the cell ID is in the list of specific IDs,
#     red if the cell is not chained (no upper or lower link),
#     and green if it is chained (has either upper or lower link).
#     """
#     if cell.cell_id in yellow_cells_ids:
#         return (1, 1, 0, 1)  # Yellow
#     elif cell.upper_link in ["None", "nan", "NaN", ""] and cell.lower_link in ["None", "nan", "NaN", ""]:
#         return (1, 0, 0, 1)  # Red
#     elif (cell.upper_link not in ["None", "nan", "NaN", ""] and
#           cell.lower_link not in ["None", "nan", "NaN", ""]) and cell.cell_id not in yellow_cells_ids:
#         return (0, 1, 0, 1)  # Green
def addAllCellsToPlot(cells, ax, ax_rng, alpha=1, show_id=False, ec='k'):
    """
        Try to optimise adding all cells to an axis
    """
    def cellPlotDict(cell):
        plot_radius = 0.5 * cell.diameter * RodShapedBacterium.sus_vis_radius_factor
        cdict = {
            'polygon': cell.getPolygon(plot_radius),
            'fc': cell.colour,
            'ec': ec,
            'alpha': alpha,
            'lw': 2 / ax_rng
        }
        return cdict

    try:
        patches = [cellPatch(**cellPlotDict(cell)) for cell in cells]
        ax.add_collection(PatchCollection(patches, match_original=True, zorder=1))
    except Exception as e:
        for cell in cells:
            cell.colour = chaining_status_to_color(cell)
            add_polygon_to_plot(ax=ax, cell=cell, fc=cell.colour, ec=ec)

def addSpringLinks(cells, ax, ax_rng, show_id=False):
    def findSprings(cell):
        for id, params in cell.springs.items():
            cell2 = cell.cell_hash[id]
            p1 = cell.rcm + params[0] * cell.ori * cell.length * 0.5
            p2 = cell2.rcm + params[1] * cell2.ori * cell2.length * 0.5
            ov = np.array(params[2:])
            p1 += cell.sus_vis_radius_factor * cell.radius * ov
            p2 -= cell2.sus_vis_radius_factor * cell2.radius * ov
            return np.array([p1[:2], p2[:2]])

    lines = [findSprings(cell) for cell in cells if cell.springs]
    anchors = [Point(pnt).buffer(0.1) for line in lines for pnt in line]
    patches = [PolygonPatch(a, fc=None, ec='k', alpha=0.8, lw=5 / ax_rng)
               for a in anchors
               ]
    connected_fraction = len([True for cell in cells if cell.springs]) / len(cells)
    print(f"{connected_fraction=}")
    ax.add_collection(
        LineCollection(lines, lw=5 / ax_rng, colors='k')
    )
    ax.add_collection(
        PatchCollection(patches, match_original=True)
    )
    return connected_fraction

def addChainLinks(cells, ax, ax_rng, show_id=False, colour='k'):
    def findChains(cell):
        try:
            upper_link = float(cell.upper_link)
            if not np.isnan(upper_link):
                cell2 = cell.cell_hash[upper_link]
                links = ChainingRodShapedBacterium.getLinkVector(cell, cell2, False)
                return [lnk[:, :-1] for lnk in links]
            else:
                print(f"NaN upper_link for cell with id: {cell.cell_id}")
        except (ValueError, KeyError) as e:
            print(f"Error processing cell with id: {cell.cell_id}, upper_link: {cell.upper_link} - {e}")
        return []

    valid_cells = [cell for cell in cells if cell.upper_link not in ["None", None, "nan", "NaN", ""] and not np.isnan(float(cell.upper_link))]

    lines = [findChains(cell) for cell in valid_cells]
    lines = [ll for line in lines if line for ll in line]

    anchors = []
    for line in lines:
        for pnt in line:
            try:
                point = Point(pnt)
                if point.is_valid:
                    anchors.append(point.buffer(0.1))
                else:
                    print("Invalid point found")
            except Exception as e:
                print("Error creating point for")

    patches = []
    for a in anchors:
        try:
            patches.append(PolygonPatch(a, fc=colour, ec=colour, alpha=0.8, lw=5 / ax_rng))
        except Exception as e:
            print("Error creating PolygonPatch")

    ax.add_collection(
        LineCollection(lines, lw=20 / ax_rng, colors=colour)
    )
    ax.add_collection(
        PatchCollection(patches, match_original=True)
    )
