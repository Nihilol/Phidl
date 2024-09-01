from phidl import Device, LayerSet, Group
from phidl import quickplot as qp # Rename "quickplot()" to the easier "qp()"
import phidl.geometry as pg
import phidl.routing as pr
import json
import numpy as np
import math
import phidl.utilities as pu
import sys
import yaml
import os
sys.path.insert(1, '..')

with open('conf.yml', 'r') as file:
    parameters = yaml.safe_load(file)

working_directory = os.getcwd()

system_params = parameters['system-params']

inner_params = parameters['inner-params']

layer_params = parameters['layer-params']

marker_params = parameters['marker-params']

dot_params = inner_params['dot-params']

barrier_params = inner_params['barrier-params']

screening_params = inner_params['screening-params']

xmax_ext_box = system_params['width_external_box']/2*1000   # Redefining so everything gets put with
xmin_ext_box = -system_params['width_external_box']/2*1000  # 0 as the center

protection_extra = system_params['protection_extra']*1000

width_bonding = system_params['width_bonding_box']*1000

transform = system_params['transform']

width_cross = barrier_params['width_cross']                      # Width of the barrier crosses

screening_gate_width = screening_params['screening_gate_width']

width_bonding_pad, length_bonding_pad = system_params['width_bonding_pad']*1000, system_params['length_bonding_pad']*1000

ymax_ext_box = 2*system_params['y_box_correction'] + system_params['height_external_box']/2*1000 # Not quite sure why the
ymin_ext_box = 2*system_params['y_box_correction'] - system_params['height_external_box']/2*1000 # the correction is needed

spacing_h_external_port = (system_params['width_external_box']*1000 - system_params['N_ports_per_edge']*system_params['width_external_port'])/(system_params['N_ports_per_edge']+1)      # Spacing between ports
spacing_v_external_port = (system_params['height_external_box']*1000 - system_params['N_ports_per_edge']*system_params['width_external_port'])/(system_params['N_ports_per_edge']+1)     # Spacing between ports

"""
Defining the layer set needed to build everything. Since this is outside the inner design
most of this is used for building the outer structure
"""

die_layer = layer_params['die_layer']
marker_layer = layer_params['marker_layer']
SiO2_layer = layer_params['SiO2_layer']
ohmic1_layer = layer_params['ohmic1_layer']
ohmic2_layer = layer_params['ohmic2_layer']
ohmic3_layer = layer_params['ohmic3_layer']
ohmic4_layer = layer_params['ohmic4_layer']
screening1_layer = layer_params['screening1_layer']
screening2_layer = layer_params['screening2_layer']
screening3_layer = layer_params['screening3_layer']
screening4_layer = layer_params['screening4_layer']
plunger1_layer = layer_params['plunger1_layer']
plunger2_layer = layer_params['plunger2_layer']
plunger3_layer = layer_params['plunger3_layer']
plunger4_layer = layer_params['plunger4_layer']
barrier1_layer = layer_params['barrier1_layer']
barrier2_layer = layer_params['barrier2_layer']
barrier3_layer = layer_params['barrier3_layer']
barrier4_layer = layer_params['barrier4_layer']
text_layer = layer_params['text_layer']
plunger_fanout_layer = layer_params['plunger_fanout_layer']
barrier_fanout_layer = layer_params['barrier_fanout_layer']
screening_fanout_layer = layer_params['screening_fanout_layer']
ohmic_fanout_layer = layer_params['ohmic_fanout_layer']
metal_buffer_layer = layer_params['metallic_buffer_layer']
ohmic_sem1_layer = layer_params['ohmic_sem1_layer']
ohmic_sem2_layer = layer_params['ohmic_sem2_layer']
ohmic_sem3_layer = layer_params['ohmic_sem3_layer']
ohmic_sem4_layer = layer_params['ohmic_sem4_layer']
screening_sem1_layer = layer_params['screening_sem1_layer']
screening_sem2_layer = layer_params['screening_sem2_layer']
screening_sem3_layer = layer_params['screening_sem3_layer']
screening_sem4_layer = layer_params['screening_sem4_layer']
barrier_sem1_layer = layer_params['barrier_sem1_layer']
barrier_sem2_layer = layer_params['barrier_sem2_layer']
barrier_sem3_layer = layer_params['barrier_sem3_layer']
barrier_sem4_layer = layer_params['barrier_sem4_layer']
plunger_sem1_layer = layer_params['plunger_sem1_layer']
plunger_sem2_layer = layer_params['plunger_sem2_layer']
plunger_sem3_layer = layer_params['plunger_sem3_layer']
plunger_sem4_layer = layer_params['plunger_sem4_layer']


lys = LayerSet()
lys.add_layer(die_layer[0], gds_layer = die_layer[1], gds_datatype = die_layer[2], alpha = 0.15)
lys.add_layer(marker_layer[0], gds_layer = marker_layer[1], gds_datatype = die_layer[2])
lys.add_layer(SiO2_layer[0], gds_layer = SiO2_layer[1], gds_datatype = SiO2_layer[2])
lys.add_layer(ohmic1_layer[0], gds_layer = ohmic1_layer[1], gds_datatype = ohmic1_layer[2])
lys.add_layer(ohmic2_layer[0], gds_layer = ohmic2_layer[1], gds_datatype = ohmic2_layer[2])
lys.add_layer(ohmic3_layer[0], gds_layer = ohmic3_layer[1], gds_datatype = ohmic3_layer[2])
lys.add_layer(ohmic4_layer[0], gds_layer = ohmic4_layer[1], gds_datatype = ohmic4_layer[2])
lys.add_layer(screening1_layer[0], gds_layer = screening1_layer[1], gds_datatype = screening1_layer[2])
lys.add_layer(screening2_layer[0], gds_layer = screening2_layer[1], gds_datatype = screening2_layer[2])
lys.add_layer(screening3_layer[0], gds_layer = screening3_layer[1], gds_datatype = screening3_layer[2])
lys.add_layer(screening4_layer[0], gds_layer = screening4_layer[1], gds_datatype = screening4_layer[2])
lys.add_layer(plunger1_layer[0], gds_layer = plunger1_layer[1], gds_datatype = plunger1_layer[2])
lys.add_layer(plunger2_layer[0], gds_layer = plunger2_layer[1], gds_datatype = plunger2_layer[2])
lys.add_layer(plunger3_layer[0], gds_layer = plunger3_layer[1], gds_datatype = plunger3_layer[2])
lys.add_layer(plunger4_layer[0], gds_layer = plunger4_layer[1], gds_datatype = plunger4_layer[2])
lys.add_layer(barrier1_layer[0], gds_layer = barrier1_layer[1], gds_datatype = barrier1_layer[2])
lys.add_layer(barrier2_layer[0], gds_layer = barrier2_layer[1], gds_datatype = barrier2_layer[2])
lys.add_layer(barrier3_layer[0], gds_layer = barrier3_layer[1], gds_datatype = barrier3_layer[2])
lys.add_layer(barrier4_layer[0], gds_layer = barrier4_layer[1], gds_datatype = barrier4_layer[2])
lys.add_layer(text_layer[0], gds_layer = text_layer[1], gds_datatype = text_layer[2])
lys.add_layer(plunger_fanout_layer[0], gds_layer = plunger_fanout_layer[1], gds_datatype = plunger_fanout_layer[2])
lys.add_layer(barrier_fanout_layer[0], gds_layer = barrier_fanout_layer[1], gds_datatype = barrier_fanout_layer[2])
lys.add_layer(screening_fanout_layer[0], gds_layer = screening_fanout_layer[1], gds_datatype = screening_fanout_layer[2])
lys.add_layer(ohmic_fanout_layer[0], gds_layer = ohmic_fanout_layer[1], gds_datatype = ohmic_fanout_layer[2])
lys.add_layer(metal_buffer_layer[0], gds_layer = metal_buffer_layer[1], gds_datatype = metal_buffer_layer[2])
lys.add_layer(ohmic_sem1_layer[0], gds_layer = ohmic_sem1_layer[1], gds_datatype = ohmic_sem1_layer[2])
lys.add_layer(ohmic_sem2_layer[0], gds_layer = ohmic_sem2_layer[1], gds_datatype = ohmic_sem2_layer[2])
lys.add_layer(ohmic_sem3_layer[0], gds_layer = ohmic_sem3_layer[1], gds_datatype = ohmic_sem3_layer[2])
lys.add_layer(ohmic_sem4_layer[0], gds_layer = ohmic_sem4_layer[1], gds_datatype = ohmic_sem4_layer[2])
lys.add_layer(screening_sem1_layer[0], gds_layer = screening_sem1_layer[1], gds_datatype = screening_sem1_layer[2])
lys.add_layer(screening_sem2_layer[0], gds_layer = screening_sem2_layer[1], gds_datatype = screening_sem2_layer[2])
lys.add_layer(screening_sem3_layer[0], gds_layer = screening_sem3_layer[1], gds_datatype = screening_sem3_layer[2])
lys.add_layer(screening_sem4_layer[0], gds_layer = screening_sem4_layer[1], gds_datatype = screening_sem4_layer[2])
lys.add_layer(barrier_sem1_layer[0], gds_layer = barrier_sem1_layer[1], gds_datatype = barrier_sem1_layer[2])
lys.add_layer(barrier_sem2_layer[0], gds_layer = barrier_sem2_layer[1], gds_datatype = barrier_sem2_layer[2])
lys.add_layer(barrier_sem3_layer[0], gds_layer = barrier_sem3_layer[1], gds_datatype = barrier_sem3_layer[2])
lys.add_layer(barrier_sem4_layer[0], gds_layer = barrier_sem4_layer[1], gds_datatype = barrier_sem4_layer[2])
lys.add_layer(plunger_sem1_layer[0], gds_layer = plunger_sem1_layer[1], gds_datatype = plunger_sem1_layer[2])
lys.add_layer(plunger_sem2_layer[0], gds_layer = plunger_sem2_layer[1], gds_datatype = plunger_sem2_layer[2])
lys.add_layer(plunger_sem3_layer[0], gds_layer = plunger_sem3_layer[1], gds_datatype = plunger_sem3_layer[2])
lys.add_layer(plunger_sem4_layer[0], gds_layer = plunger_sem4_layer[1], gds_datatype = plunger_sem4_layer[2])



layer_marker =  lys['marker']
layer_protection_bonding_pad = lys['SiO2']
layer_die = lys['die']
layer_text = lys['text']

layer_dummy_bonds = [lys['Ohmic fanout'], lys['Barrier fanout'], lys['Plunger fanout']]
layers_outer = [lys['Plunger fanout'], lys['Barrier fanout'], lys['Screening fanout'], lys['Ohmic fanout']]
layer_dot_array = [lys['Plunger Ti/Pd 1'], lys['Plunger Ti/Pd 2'], lys['Plunger Ti/Pd 3'], lys['Plunger Ti/Pd 4']]
layer_screening_array = [lys['Screening Ti/Pd 1'], lys['Screening Ti/Pd 2'], lys['Screening Ti/Pd 3'], lys['Screening Ti/Pd 4']]
layer_barrier_array = [lys['Barrier Ti/Pd 1'], lys['Barrier Ti/Pd 2'], lys['Barrier Ti/Pd 3'], lys['Barrier Ti/Pd 4']]
layer_ohmic_array = [lys['Ohmic Al 1'], lys['Ohmic Al 2'], lys['Ohmic Al 3'], lys['Ohmic Al 4']]

layer_dot_SEM_array = [lys['Plunger Ti/Pd SEM1'], lys['Plunger Ti/Pd SEM2'], lys['Plunger Ti/Pd SEM3'], lys['Plunger Ti/Pd SEM4']]
layer_screening_SEM_array = [lys['Screening Ti/Pd SEM1'], lys['Screening Ti/Pd SEM2'], lys['Screening Ti/Pd SEM3'], lys['Screening Ti/Pd SEM4']]
layer_barrier_SEM_array = [lys['Barrier Ti/Pd SEM1'], lys['Barrier Ti/Pd SEM2'], lys['Barrier Ti/Pd SEM3'], lys['Barrier Ti/Pd SEM4']]
layer_ohmic_SEM_array = [lys['Ohmic Al SEM1'], lys['Ohmic Al SEM2'], lys['Ohmic Al SEM3'], lys['Ohmic Al SEM4']]

def write_SCON(name, coord_list, size):

    if not os.path.isdir('Marker SCONS'):
        os.mkdir('Marker SCONS')

    file_name = "Marker SCONS/" + name + ".SCON"
    with open(file_name, "w") as f:
        first_line = "/* xxxxxxxxxxxxxxxxxxxx " + name + ".SCON xxxxxxxxxxxxxxxxxxxx */\n"
        f.write(first_line)
        f.write("VER,1.120615\n")
        f.write("MODEL,ELS-F125-50B\n")
        f.write("ID,F125PRU\n")
        f.write("MF,1000\n")
        f.write("FLD,D:\\\n")
        f.write("C1,1,1,0.000000\n")
        f.write("C2,1,1,0.000000\n")
        f.write("C3,1,1,0.000000\n")
        f.write("C4,1,1,0.000000\n")
        f.write("C5,1,1,0.000000\n")
        f.write("C6,1,1,0.000000\n")
        f.write("C7,1,1,0.000000\n")
        f.write("!START\n")
        f.write("PX,200000\n")
        f.write("SZ,500\n")
        f.write("STYP,1\n")
        if size == "local":
            coord_string = "RL4,"
        if size == "global":
            coord_string = "RS,"
        for i in coord_list:
            formatted_x = "%.10f" % i[0]
            formatted_y = "%.10f" % i[1]
            if i == coord_list[-1]:
                coord_string = coord_string + str(formatted_x) + "," + str(formatted_y)
            else:
                coord_string = coord_string + str(formatted_x) + "," + str(formatted_y) + ","
        coord_string = coord_string + "\n"
        f.write(coord_string)
        f.write("!END\n")
        f.write("EST1,12.0000\n")
        f.write("EST2,24.0000\n")
        f.write("EST3,60.4500\n")
        f.write("EST4,30.2200\n")
        f.write("EST5,0.4000\n")
        f.write("EST6,12")


def place_sem_structure(device, marker_structure, cardinality,
                        xmin, ymin):

    """
    place_sem_structure is a function used to place the smaller structures around
    the chip, that we can use as dummy structures to perform SEM.

    :device: refers to the device the sem structure is added to

    :marker_structure: is the inner structure, currently without inner markers added.

    :cardinality: is a deprecated keyword, which is still necesarry to place them,
    but it no longer serves it original purpose. It used to be that this function required
    a reference device, from which you could then place the SEM structure via cardinality.

    :xmin: and ymin are simple coordinates used to place the structure.
    """

    sw_marker = inner_marker_generator(2*system_params['width_inner_marker'], 2*system_params['length_inner_marker']*1000, layer_marker, cardinal = "A")
    se_marker = inner_marker_generator(2*system_params['width_inner_marker'], 2*system_params['length_inner_marker']*1000, layer_marker, cardinal = "B")
    nw_marker = inner_marker_generator(2*system_params['width_inner_marker'], 2*system_params['length_inner_marker']*1000, layer_marker, cardinal = "C")
    ne_marker = inner_marker_generator(2*system_params['width_inner_marker'], 2*system_params['length_inner_marker']*1000, layer_marker, cardinal = "D")


    if cardinality == "SE":
        ref_dummy_devices = device.add_ref(marker_structure)
        ref_dummy_devices.xmin = xmin
        ref_dummy_devices.ymin = ymin

    if cardinality == "SW":
        ref_dummy_devices = device.add_ref(marker_structure)
        ref_dummy_devices.ymin = ymin
        ref_dummy_devices.xmin = xmin

    if cardinality == "NE":
        ref_dummy_devices = device.add_ref(marker_structure)
        ref_dummy_devices.ymin = ymin
        ref_dummy_devices.xmin = xmin

    if cardinality == "NW":
        ref_dummy_devices = device.add_ref(marker_structure)
        ref_dummy_devices.ymin = ymin
        ref_dummy_devices.xmin = xmin
    
    sw_marker.move(destination = (ref_dummy_devices.center[0] - 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] - 0.5*ref_dummy_devices.ysize))
    se_marker.move(destination = (ref_dummy_devices.center[0] + 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] - 0.5*ref_dummy_devices.ysize))
    nw_marker.move(destination = (ref_dummy_devices.center[0] - 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] + 0.5*ref_dummy_devices.ysize))
    ne_marker.move(destination = (ref_dummy_devices.center[0] + 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] + 0.5*ref_dummy_devices.ysize))
    
    device.add_ref(sw_marker)#.move(destination = (ref_dummy_devices.center[0] - 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] - 0.5*ref_dummy_devices.ysize))
    device.add_ref(se_marker)#.move(destination = (ref_dummy_devices.center[0] + 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] - 0.5*ref_dummy_devices.ysize))
    device.add_ref(nw_marker)#.move(destination = (ref_dummy_devices.center[0] - 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] + 0.5*ref_dummy_devices.ysize))
    device.add_ref(ne_marker)#.move(destination = (ref_dummy_devices.center[0] + 0.5*ref_dummy_devices.xsize, ref_dummy_devices.center[1] + 0.5*ref_dummy_devices.ysize))


    return sw_marker, se_marker, nw_marker, ne_marker

def create_and_place_dummy_bond_pads(device, amount_dummy_bonding_pads, width, length, spacing, direction, destination, layer_dummy_bonds):
    
    """ 
    The function below is used to generate the dummy bond pads, on the outer skirts
    of the total design.

    :device: refers to the device the dummy pads will be added to.

    :amount_dummy_bonding_pads: refers to how many dummy pads there will be generated.

    :width: is the width of the dummy bond pads.

    :length: is the length of the dummy bond pad.

    :spacing: is the spacing in between each dummy bond pad.

    :direction: is whether the grouping will be directed horizontally or vertically.

    :destination: is the coordinate the grouping will be moved to
    """
    
    buffer_params = parameters["buffer-params"]
    
    buffer_radius = buffer_params["circle_radius"]*1000
    
    dummy_bonding_pads = []
    buffer_pads_locs = []
    j = 0
    for i in range(amount_dummy_bonding_pads):
        if j == 3:
            j = 0
        if j > 0:
            if direction == 'y':
                buffer_pads_locs.append([destination[0] + 0.5*width, destination[1] + 0.5*length + 2*i*spacing])
            else:
                buffer_pads_locs.append([destination[0] + 0.5*width + 2*i*spacing, destination[1] + 0.5*length])
        dummy_bonding_pads.append(device << pg.rectangle(size = (length, width), layer = layer_dummy_bonds[j]))
        j += 1
    
    dummy_bonding_pads_group = Group(dummy_bonding_pads)
    dummy_bonding_pads_group.distribute(direction = direction, spacing = spacing)
    dummy_bonding_pads_group.move(destination = destination)

    buffer_pads = []    

    for i in range(len(buffer_pads_locs)):
        buffer_pads.append("buffer_pad" + str(i))
        buffer_pads[i] = device << pg.circle(radius = buffer_radius, layer = metal_buffer_layer[1])
        buffer_pads[i].center = dummy_bonding_pads[i + 1].center

    return dummy_bonding_pads


def inner_marker_generator(width, length, layer, cardinal = ""):

    """
    :inner_marker_generator: is used to generate the inner, smaller markers, both used
    on the inside of the device, but also inside the larger markers on the outskirts
    of the design.

    :width: is the width of each cross.

    :length: is the length of each cross. 

    :cardinal: is used as an identifying marker. It can be empty, and will then not
    have any effect, however, it is currently used as an alpha-numerical ordering code,
    such that A, B, C, D are the lower left, lower right, upper left, and upper right
    corners of the device, respectively. The numbers are assigned such that each corner
    cross signifies the 0 position, and all other crosses in that alpha section are
    given relative coordinates to that.
    """

    marker = Device()

    cross_inner_marker = marker << pg.cross(length = length, width = width, layer = layer)
    
    if cardinal:
        cardinal_text = marker << pg.text(cardinal, layer = layer, size = 400)
        cardinal_text.move(destination = (-1500, 1200))   

    inner_marker_joined = pg.union(marker, by_layer = True)

    return inner_marker_joined


def global_marker_generator():
    
    """ 
    Generate the outer markers, which are used on the outskirt of the design, for rough
    alignment. 
    No parameters are given; if wanting to change their size, look at the source-code.
    The pg.boolean function is used to cut out a specfic size inside the outer marker,
    such that an inner marker can be placed there.
    """
    
    cross_outer_marker =  pg.cross(length = system_params['length_outer_marker']*1000, width = system_params['width_outer_marker']*1000, layer = layer_marker)

    outer_marker_device_joined = pg.union(cross_outer_marker, by_layer = True)

    rect_box_cut_outer_marker = pg.rectangle(size = (system_params['edge_box_cut_outer_marker']*1000, system_params['edge_box_cut_outer_marker']*1000), layer = layer_marker)
    rect_box_cut_outer_marker.center = [0,0]

    outer_marker_joined_cut = pg.boolean(outer_marker_device_joined,rect_box_cut_outer_marker,'A-B', layer = layer_marker)

    return outer_marker_joined_cut


def marker_device_generator(device, x, y, cardinal):
    
    """ 
    This function combines the functionality of the inner_marker_generator and the
    global_marker_generator, to properly generate the inner/outer markers used on the
    outskirts of the design.

    :x: is the x direction numerical ordering code mentioned above. It is given relatively
    to the corner marker on each side of the device. 

    :y: is the y direction numerical ordering code.

    :cardinal: is the absolute placement on the chip. Not functionally needed to place 
    the chips, but needed for the text-generation.
        A = South west
        B = South east
        C = North west
        D = North east
    """

    layer = lys['marker']

    amount_of_markers = len(x)
    marker_spacing = 200000

    marker_device = Device()

    inner_markers, outer_markers = [], []
    position = []
    cardinal_list = []
    names = []

    for i in range(amount_of_markers):
        position.append(str(x[i]) + str(y[i]))
        position[i] = marker_device << pg.text((str(x[i]) + str(y[i])), layer = layer, size = 400)
        cardinal_list.append(cardinal)
        cardinal_list[i] = marker_device << pg.text(cardinal, layer = layer, size = 400)
        names.append(cardinal + str(x[i]) + str(y[i]))
        inner_markers.append("inner_left_" + str(i))
        outer_markers.append("outer_left_" + str(i))
        inner_markers[i] = inner_marker_generator(width = system_params['width_inner_marker'], length = system_params['length_inner_marker']*1000, layer = layer)
        outer_markers[i] = marker_device << global_marker_generator()

    def marker_mover(list1, list2, i, x, y):
        list1[i].move(destination = (x*marker_spacing, y *marker_spacing))
        list2[i].move(destination = (x*marker_spacing, y *marker_spacing))

    for j in range(amount_of_markers):
        marker_mover(inner_markers, outer_markers, j, x[j], y[j])
        position[j].move(destination = (x[j]*marker_spacing + 700, y[j]*marker_spacing + 700))
        cardinal_list[j].move(destination = (x[j]*marker_spacing - 1200, y[j]*marker_spacing + 700))

    device.add_ref(inner_markers)


    inner_markers_group = Group(inner_markers)

    return marker_device, inner_markers, inner_markers_group, names

""" 
Credit to Fabrizio Berritta for the two following functions.

They calculate the midpoint for each of the bonding pads, simply from
giving them the integer number of placements which the bond-pad is placed on.
"""

def x_midpoint_ext_port_h(N):
    x0 = spacing_h_external_port*N+system_params['width_external_port']*(N-1)+ system_params['width_external_port']/2 + xmin_ext_box
    return x0
def y_midpoint_ext_port_v(N):
    y0 = spacing_v_external_port*N+system_params['width_external_port']*(N-1)+ system_params['width_external_port']/2 + ymin_ext_box
    return y0

def calculate_placements(amount_of_dots):
        
        """ 
        The function calculates the amount of routings needed to the south
        and north side based on the amount of dots input.
        It is functional for everything below 10 dots, but
        after that, more functionality is needed.
        It is based on an alternating counter, since each inner
        cross barrier is alternated between the north and south
        routing.
        
        :amount_of_dots: The amount of dots in the design.
        """
        
        north_counter = 0
        south_counter = 0
        if amount_of_dots > 3:
            for i in range(3, amount_of_dots + 1):
                if i%2 == 0:
                    north_counter += 1
                    south_counter += 1
        if amount_of_dots == 1:
            north_addition = 2
            south_addition = 2
        elif amount_of_dots%2 == 0:
            north_addition = 2
            south_addition = 2
        elif amount_of_dots%2 != 0:
            south_addition = 2
            north_addition = 3
        north_end = 2*amount_of_dots + north_addition + north_counter
        south_end = 2*amount_of_dots + south_addition + south_counter
        west_end = 6
        east_end = 6
        
        north_left = 12 - north_end
        south_left = 12 - south_end
        west_left = 6
        east_left = 6
        
        north_placement = np.arange(north_left/2 + 1, north_end + north_left/2 + 1, dtype = int)
        south_placement = np.arange(south_left/2 + 1, south_end + south_left/2 + 1, dtype = int)
        west_placement = np.arange(west_left/2 + 1, west_end + west_left/2 + 1, dtype = int)
        east_placement = np.arange(east_left/2 + 1, east_end + east_left/2 + 1, dtype = int)
        
        return north_placement, south_placement, west_placement, east_placement

def make_appendage(object, cardinality, appendage_lengths,  layer_barrier, layer, port_routing_direction = [0, 0]):
        
        """ 
        A function to generate the appendage gates that are on the cross barriers,
        so that they can be routed out properly, and not directly from where they are
        in the center of the design. 
        
        :object: The object that the appendage is attached to.
        
        :cardinality: The cardinality of the appendage. 1 is down, 3 is up.
        
        :appendage_lengths: The length of the appendage.
        
        :port_routing_direction: The skew of each appendage port, so ensure no
        thinning of gates.
        """
        
        if cardinality == '3':
            appendage = pg.polygon_ports(xpts = [object[0].center[0] - 0.5*width_cross,
                                                 object[0].center[0] + 0.5*width_cross, 
                                                 object[0].center[0] + 0.5*width_cross, 
                                                 object[0].center[0] - 0.5*width_cross], 
             ypts = [                            object[0].ymax, 
                                                 object[0].ymax, 
                                                 object[0].ymax + appendage_lengths + port_routing_direction[0],
                                                 object[0].ymax + appendage_lengths + port_routing_direction[1]], layer = layer_barrier[layer])
            appendage_port = appendage.ports['3']
        else:
            appendage = pg.polygon_ports(xpts = [object[0].center[0] - 0.5*width_cross,
                                                 object[0].center[0] + 0.5*width_cross, 
                                                 object[0].center[0] + 0.5*width_cross, 
                                                 object[0].center[0] - 0.5*width_cross], 
             ypts = [                            object[0].ymin - appendage_lengths - port_routing_direction[0], 
                                                 object[0].ymin - appendage_lengths - port_routing_direction[1], 
                                                 object[0].ymin,
                                                 object[0].ymin], layer = layer_barrier[layer])
            appendage_port = appendage.ports['1']
        return appendage, appendage_port


def protection_pads(device, width, length, destination):
   
    """ 
    :protection_pads: is a function which builds the protection pads underneath
    the other bonding pads, except for the ohmics (user-defined).

    :device: refers to the device they are added onto.

    :width: is the width of the protection pad.

    :length: is the length of the protection pad.

    :destination: is the coordinate on the chip where the protection pad will be placed.
    """
    
    protection_pad = device << pg.rectangle(size = (width, length), layer = layer_protection_bonding_pad)
    protection_pad.move(destination = destination)
    return protection_pad


def cross_barriers(mid_point, width, length, layer, port_cardinal, port_routing_direction):
    
    """ 
    :cross_barriers: is the function that builds the cross barriers seen on the inside
    of the 2xN device. 

    :mid_point: is a touple with x and y coordinates.

    :width: is the width of each cross.

    :length: is the length of each cross.

    :layer: is the layer each cross is assigned to.

    :port_cardinal: is the designation for whether you the output port to be on the
    north or south side of the cross.

    :port_routing_direction: is a list of length two, with numerical values, which
    is meant to add a skew to the chosen output port. This is done so that no
    thinning of the gates occur, as they fan into the inner structure.
    """

    if port_cardinal == '3':
        vertical = pg.polygon_ports(xpts = [        mid_point[0] - 0.5*width,
                                                    mid_point[0] + 0.5*width, 
                                                    mid_point[0] + 0.5*width, 
                                                    mid_point[0] - 0.5*width], 
                               ypts = [             mid_point[1] - 0.5*length, 
                                                    mid_point[1] - 0.5*length, 
                                                    mid_point[1] + 0.5*length + port_routing_direction[1],
                                                    mid_point[1] + 0.5*length + port_routing_direction[0]], layer = layer)
        horizontal = pg.polygon_ports(xpts = [      mid_point[0] - 0.5*length,
                                                    mid_point[0] + 0.5*length, 
                                                    mid_point[0] + 0.5*length, 
                                                    mid_point[0] - 0.5*length], 
                               ypts = [             mid_point[1] - 0.5*width, 
                                                    mid_point[1] - 0.5*width, 
                                                    mid_point[1] + 0.5*width,
                                                    mid_point[1] + 0.5*width], layer = layer)
        port = vertical.ports[port_cardinal]
    elif  port_cardinal == '1':
        vertical = pg.polygon_ports(xpts = [        mid_point[0] - 0.5*width,
                                                    mid_point[0] + 0.5*width, 
                                                    mid_point[0] + 0.5*width, 
                                                    mid_point[0] - 0.5*width], 
                               ypts = [             mid_point[1] - 0.5*length - port_routing_direction[1], 
                                                    mid_point[1] - 0.5*length - port_routing_direction[0], 
                                                    mid_point[1] + 0.5*length,
                                                    mid_point[1] + 0.5*length], layer = layer)
        horizontal = pg.polygon_ports(xpts = [      mid_point[0] - 0.5*length,
                                                    mid_point[0] + 0.5*length, 
                                                    mid_point[0] + 0.5*length, 
                                                    mid_point[0] - 0.5*length], 
                               ypts = [             mid_point[1] - 0.5*width, 
                                                    mid_point[1] - 0.5*width, 
                                                    mid_point[1] + 0.5*width,
                                                    mid_point[1] + 0.5*width], layer = layer)

        port = vertical.ports[port_cardinal]
    else:
        raise Exception('The port cardinal direction you entered was not valid. Please make it 1 for the south port, and 3 for the north port.')
    return (horizontal, vertical), port

def barrier_creation(dot, dot_radius, width_barrier, length_barrier, layer, rotations, routing_direction):
    """ 
    :barrier_creation: Creates the barriers that are in between the sensor dots
    and the regular dots. It is only used in relation to the creation of the sensor dots
    and the positioning of this barrier is only relative to the sensor dots.

    :dot: is the dot from which the relative coordinates will be derived.

    :dot_radius: Is the radius of the dot above, which is used to determine the relative
    coordinates for the barrier.

    :width_barrier: is the width of the barrier.

    :length_barrier: is the length of the barrier.

    :layer: is the layer the barrier is assigned to.

    :rotations: is the rotation of the barrier, relative to the dot. For this specific 
    purpose, all rectangle barriers are created such that they are aligned vertically,
    and then rotated around the dot, such that they are aligned horizontally.

    :routing_direction: is either "left" or "right", meaning which sensor dot the barrier
    belongs to.
    """

    distance_from_sensordot_to_barrier = dot_params['distance_from_sensordot_to_barrier']

    if routing_direction == 'left':
        if rotations == 90:
            polygon = pg.polygon_ports(xpts = [             dot.center[0] + dot_radius + distance_from_sensordot_to_barrier,
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier], 
                                    ypts = [                dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] + 0.25*length_barrier,
                                                            dot.center[1] + 0.25*length_barrier + 0.15*length_barrier], layer = layer)
        if rotations == -90:
            polygon = pg.polygon_ports(xpts = [             dot.center[0] + dot_radius + distance_from_sensordot_to_barrier,
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier], 
                                    ypts = [                dot.center[1] - 0.25*length_barrier - 0.15*length_barrier, 
                                                            dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] + 0.25*length_barrier,
                                                            dot.center[1] + 0.25*length_barrier], layer = layer)
    else:
        if rotations == 90:
            polygon = pg.polygon_ports(xpts = [             dot.center[0] + dot_radius + distance_from_sensordot_to_barrier,
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier], 
                                    ypts = [                dot.center[1] - 0.25*length_barrier - 0.15*length_barrier, 
                                                            dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] + 0.25*length_barrier,
                                                            dot.center[1] + 0.25*length_barrier], layer = layer)
        if rotations == -90:
            polygon = pg.polygon_ports(xpts = [             dot.center[0] + dot_radius + distance_from_sensordot_to_barrier,
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier + width_barrier, 
                                                            dot.center[0] + dot_radius + distance_from_sensordot_to_barrier], 
                                    ypts = [                dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] - 0.25*length_barrier, 
                                                            dot.center[1] + 0.25*length_barrier,
                                                            dot.center[1] + 0.25*length_barrier + 0.15*length_barrier], layer = layer)
    return polygon, polygon.ports


def dot_creation(gate_additional_length, center, radius, angle_resolution, layer, width_gate,
                left_tilt, right_tilt, angle_of_gate, barrier_layer,
                 ports, port_routings, barrier_rot = [], mirror = False):
    
    """ 
    :dot_creation: Creates the dots that are used in the dot array.

    :gate_additional_length: is the additional length of the gate, which is used
    to ensure proper fan-out, even as the amount of dots become high.

    :center: is the coordinate for the center of the dot.

    :radius: is the radius of the dot.

    :angle_resolution: is the resolution of the angles in the gds file.

    :layer: is the layer the dot is assigned to.

    :width_gate: is the width of the gate.

    :left_tilt: is the assigned left skew on the gate, meant to ensure that no thinning
    of the gate happens as it is fanned out.

    :right_tilt: is the assigned right skew on the gate, meant to ensure that no thinning
    of the gate happens as it is fanned out.

    :angle_of_gate: is the angle by which the gate will be rotated around the dot.
    This is needed, since the gates need to fan out to the north and the south.

    :barrier_layer: This is currently only used for the creation of the sensor dots.
    It is the layer that the barriers created via the above function will be assigned to.

    :barrier_port_direction: This is currently only used for the creation of
    the sensor dots. It is the direction that the barriers will be fanned out of.

    :barrier_rot: This is currently only used for the creation of the sensor dots. 
    It is the rotation of the barriers, relative to the sensor dot.

    :routing_directions: This is currently only used for the creation of the sensor dots.
    It is either "left" or "right" meaning which sensor dot the barrier belongs to.

    :mirror: This is used to mirror the gates. This is so that each dot gate can be
    made in the same manner, and then as they are rotated around, to be fanned out in the
    opposite direction, they can be mirrored to ensure that the gates fan out in the
    correct fashion.
    """
    
    dot = pg.circle(radius = radius, angle_resolution = angle_resolution, layer = layer)
    dot.center = center

    gate = pg.polygon_ports(xpts = [dot.center[0] - 0.5*width_gate,
                                                    dot.center[0] + 0.5*width_gate, 
                                                    dot.center[0] + 0.5*width_gate, 
                                                    dot.center[0] - 0.5*width_gate], 
                                            ypts = [dot.center[1] + 0.90*radius, 
                                                     dot.center[1] + 0.90*radius, 
                                                     dot.center[1] + 0.95*radius + gate_additional_length + right_tilt,
                                                     dot.center[1] + 0.95*radius + gate_additional_length + left_tilt], layer = layer)
    gate.rotate(angle = angle_of_gate, center = dot.center)
    if mirror:
        gate.mirror(p1 = gate.center, p2 = (gate.center[0], gate.center[1] + 1))

    barrier_params = inner_params['barrier-params']

    width_barrier = barrier_params['width_barrier']
    length_barrier = 200

    barrier_array = []
    
    barrier_ports = []

    j = 0
    
    gate_port = '3'

    if len(barrier_rot) != 0:
        gate.rotate(angle = 180, center = gate.center)
        gate.mirror(p1 = [gate.xmin, gate.center[1]], p2 = [gate.xmax, gate.center[1]])
        for i in barrier_rot:
            # barrier_i, full_barrier_ports = barrier_creation(dot, radius, width_barrier, length_barrier, layer = barrier_layer, rotations = i,
            #                                                  routing_direction = routing_directions)
            barrier_i, full_barrier_ports = cross_barriers(mid_point = [dot.center[0] + radius + 100/4, dot.center[1]],
                                                                    width = 45,
                                                                    length = 82, layer = barrier_layer,
                                                                    port_cardinal = ports[j],
                                                                    port_routing_direction = port_routings[j])
            barrier_array.append(barrier_i)
            barrier_array[j][0].rotate(angle = i, center = dot.center)
            barrier_array[j][1].rotate(angle = i, center = dot.center)
            barrier_ports.append(full_barrier_ports)
            j += 1
    return dot, gate, barrier_array, gate.ports[gate_port], barrier_ports


def build_2xN(layer, amount_of_dots, layer_dot, layer_barrier, layer_screening, layer_ohmic,
              layers_for_outer, local_markers_pos, dose_test = False):
    
    """
    :inner_params: is a dictionary which contains all the parameters for the design. 
    It is read from the yaml file conf.yml
    
    :layer: is a parameter meant to make it easier to create more of these design onto a single
    chip. For this purpose, if you want to create 4 of these designs onto a single chip,
    you create a 4-long list of layers for each subgroup of layers. This is done 
    so that in fabrication, the layers are not confused with each other. This parameter
    is an integer representing the 1st, 2nd, 3rd, or 4th design on the chip.
    
    :amount_of_dots_ is the amount of dots in a single row. So if you give it 3, it will
    create a 2x3 design.
    
    :layer_dot: is a 4-long list of the layers meant for the dots.
    
    :layer_barrier: is a 4-long list of the layers meant for the barriers.
    
    :layer_screening: is a 4-long list of the layers meant for the screening gates.
    
    :layer_ohmic: is a 4-long list of the layers meant for the ohmic contacts.
    
    :layers_for_outer: is a 4-long list of the layers meant for the outermost gates.
    
    :name_of_chip: is the name of the chip, which is used to name the gds file.
    
    :dose_test: is a boolean which is used to determine whether it saves a file of the
    inner device without markers, for dose testing purposes. False by default.
    """
    
    # Define all constants and parameters from params dictionary ------------------------------
    
    dot_params = inner_params['dot-params']
    screen_params = inner_params['screening-params']
    barrier_params = inner_params['barrier-params']
    
    dot_radius = dot_params['dot_radius']                        # Dot radius
    sensordot_radius = dot_params['sensordot_radius']            # Sensordot radius
    dot_width_gate = dot_params['dot_width_gate']                # width of connection dot gate
    angle_resolution = dot_params['angle_resolution']            # Resolution of design
    sensor_width_gate = dot_params['sensor_width_gate']          # Width of connecting sensor dot gate
    spacing_dots = dot_params['pitch_dots']                    # Spacing between dots
    start_x = dot_params['start_x']                              # The starting x-coordinate dot 1,1
    start_y = dot_params['start_y']                              # Starting y-coordinate of dot 1,1
    screen_length = screen_params['screen_length']                  # Size of the screening gate above dots
    gate_length_screen = screen_params['gate_length_screen']        # Length of the connecting screening gate
    width_barrier = barrier_params['width_barrier']                  # Width of the barriers
    length_cross = barrier_params['length_cross']                    # Length of the barrier crosses
    middle_cross_length_addition = barrier_params['middle_cross_length_addition'] # Same as above, but for the middle barrier crosses
    gate_length_start = inner_params['gate_length_start']          # Starting point for list of increasingly large gate lengths, so ensure proper fan-out
    gate_length_end = inner_params['gate_length_end']              # End point for list of increasingly large gate lengths, so ensure proper fan-out
    middle_gate_addition = inner_params['middle_gate_addition']    # If uneven number of gates, the middle gate will be slightly larger
    appendage_length_addition = inner_params['appendage_length_addition']    # Same as above, but for the appendages
    overhang_screening = screen_params['overhang_from_plunger_gates']      # Overhang of the screening gate
    length_dot_perimeter_to_screening = screen_params['length_from_dot_perimeter']
    distance_from_sensordot_to_screening_gate = dot_params['distance_from_sensordot_to_screening_gate']
    distance_from_sensordot_to_dot = dot_params['distance_from_sensor_to_dot']
    lenght_barrier_between_sensor_and_dot = barrier_params['lenght_barrier_between_sensor_and_dot']
    tilt_barrier_between_sensor_and_dot = barrier_params['tilt_barrier_between_sensor_and_dot']
    
    
    # ----------------------------------------------------------------------------------------
    
    """
    The following code creates empty arrays, which are filled with the values
    of the skew meant to be applied to the gates. The skew is meant to compensate
    for the thinning of gates as they are fanned out at non-zero degree angles.
    """
    
    
    dot_tilts_right_2xN = np.zeros((2, amount_of_dots), dtype = float)      
    dot_tilts_left_2xN = np.zeros((2, amount_of_dots), dtype = float)
    dot_rotations_2xN = np.zeros((2, amount_of_dots), dtype = float)
    
    tilt_start = -50
    tilt_end = 20
    
    dot_tilts_left_2xN[0] = np.linspace(tilt_start, tilt_end, num = amount_of_dots)
    dot_tilts_left_2xN[1] = np.linspace(tilt_start, tilt_end, num = amount_of_dots)
    dot_tilts_right_2xN[0] = np.linspace(tilt_end, tilt_start, num = amount_of_dots)
    dot_tilts_right_2xN[1] = np.linspace(tilt_end, tilt_start, num = amount_of_dots)
    
    if amount_of_dots == 3:
        dot_tilts_left_2xN[0][1] = -50
        dot_tilts_right_2xN[0][1] = 0
        dot_tilts_left_2xN[1][1] = -50
        dot_tilts_right_2xN[1][1] = 0
    
    dot_rotations_2xN[0] = np.linspace(180, 180, num = amount_of_dots)
    dot_rotations_2xN[1] = np.linspace(0, 0, num = amount_of_dots)
    
    sensor_tilts_right_2x2 = [30, 0]

    sensor_tilts_left_2x2 = [0, 30]

    sensor_rotations_2x2 = [90, -90]

    """
    The following code initialises the additional lengths, added to the various
    gates, such that proper fan-out is achieved even at higher number of dots. 
    """    

    middle_gate_length = gate_length_end + middle_gate_addition
    
    appendage_length_start = gate_length_start + appendage_length_addition
    appendage_length_end = gate_length_end + appendage_length_addition
    appendage_middle_gate_length = middle_gate_length + appendage_length_addition

    middle_cross_appendage_length_start = appendage_length_start
    middle_cross_middle_gate_length = appendage_middle_gate_length 
    middle_cross_appendage_length_end = appendage_length_end + middle_cross_length_addition
    
    if (amount_of_dots - 1)%2 == 0:
        appendage_length_addition_1 = np.linspace(appendage_length_start, appendage_length_end, num = round(np.floor((amount_of_dots - 1)/2)))
        appendage_length_addition_2 = np.linspace(appendage_length_end, appendage_length_start, num = round(np.floor((amount_of_dots - 1)/2)))
        appendage_length_addition = np.hstack((appendage_length_addition_1, appendage_length_addition_2))
    else:
        appendage_length_addition_1 = np.linspace(appendage_length_start, appendage_length_end, num = round(np.floor((amount_of_dots - 1)/2)))
        appendage_length_addition_2 = np.linspace(appendage_middle_gate_length, appendage_middle_gate_length, 1)
        appendage_length_addition_3 = np.linspace(appendage_length_start, appendage_length_end, num = round(np.floor((amount_of_dots - 1)/2)))[::-1]
        appendage_length_addition = np.hstack((appendage_length_addition_1, appendage_length_addition_2, appendage_length_addition_3))
    
    
    if amount_of_dots%2 == 0:
        gate_length_addition_1 = np.linspace(gate_length_start, gate_length_end, num = round(np.floor(amount_of_dots/2)))
        gate_length_addition_2 = np.linspace(gate_length_start, gate_length_end, num = round(np.floor(amount_of_dots/2)))[::-1]
        gate_length_addition = np.hstack((gate_length_addition_1, gate_length_addition_2))
        appendage_length_middle_addition_1 = np.linspace(middle_cross_appendage_length_start, middle_cross_appendage_length_end, num = round(np.floor((amount_of_dots)/2)))
        appendage_length_middle_addition_2 = np.linspace(middle_cross_appendage_length_end, middle_cross_appendage_length_start, num = round(np.floor((amount_of_dots)/2)))
        appendage_length_middle_addition = np.hstack((appendage_length_middle_addition_1, appendage_length_middle_addition_2))
    else:
        gate_length_addition_1 = np.linspace(gate_length_start, gate_length_end, num = round(np.floor(amount_of_dots/2)))
        gate_length_addition_2 = np.linspace(middle_gate_length, middle_gate_length, 1)
        gate_length_addition_3 = np.linspace(gate_length_start, gate_length_end, num = round(np.floor(amount_of_dots/2)))[::-1]
        gate_length_addition = np.hstack((gate_length_addition_1, gate_length_addition_2, gate_length_addition_3))
        appendage_length_middle_addition_1 = np.linspace(middle_cross_appendage_length_start, middle_cross_appendage_length_end, num = round(np.floor((amount_of_dots)/2)))
        appendage_length_middle_addition_2 = np.linspace(middle_cross_middle_gate_length, middle_cross_middle_gate_length, 1)
        appendage_length_middle_addition_3 = np.linspace(middle_cross_appendage_length_start, middle_cross_appendage_length_end, num = round(np.floor((amount_of_dots)/2)))[::-1]
        appendage_length_middle_addition = np.hstack((appendage_length_middle_addition_1, appendage_length_middle_addition_2, appendage_length_middle_addition_3))
    
    sensor_dot_dot_radii_difference = sensordot_radius - dot_radius # Difference between the two sizes, so that relative spacing can be achieved

    """
    Define dots and their positions, including initialising arrays for their ports and positions
    """
    dot_positions_2xN = np.zeros((2, amount_of_dots), dtype = tuple)
    dot_2xN = np.zeros((2, amount_of_dots), dtype = tuple)
    dot_ports_2xN = np.zeros((2, amount_of_dots), dtype = tuple)
    
    for i in range(amount_of_dots):
        dot_positions_2xN[0, i] = [start_x + i*spacing_dots, start_y]
        dot_positions_2xN[1, i] = [start_x + i*spacing_dots, start_y + spacing_dots]

    
    """
    Define the sensor dots and their positions, including initialising arrays for their ports and positions
    """
    
    sensordot_positions_2xN = [[start_x - distance_from_sensordot_to_dot, start_y],
                               [start_x + (amount_of_dots-1)*spacing_dots + distance_from_sensordot_to_dot,
                                start_y + spacing_dots]]    
    sensor_dots_2xN = np.zeros((2), dtype = tuple)
    
    sensor_dot_ports_2xN = np.zeros((2), dtype = tuple)
    
    barrier_rotations_2x2 = [-90, 90]
    
    sensor_routing_directions = ['left', 'right']
    
    barrier_port_direction = [('1', '3'), ('3', '1')] # The correct ports for the barriers on the side of the sensordots
    
    embedded_barrier_ports_2x2 = []
    
    dot_gates = np.zeros((2, amount_of_dots), dtype = tuple)
    
    """
    Build the actual dots based on the positions and parameters defined above.
    """
    
    D = Device()
    
    for i in range(amount_of_dots):
        dot_2xN[0, i], dot_gates[0, i], barrier_array_1, dot_ports_2xN[0, i], _ = dot_creation(gate_length_addition[i], 
                            dot_positions_2xN[0, i], dot_radius, angle_resolution,
                            layer_dot[layer], dot_width_gate, dot_tilts_left_2xN[0, i],
                            dot_tilts_right_2xN[0, i], dot_rotations_2xN[0, i], layer_barrier[layer], [], [], mirror = True)
        dot_2xN[1, i], dot_gates[1, i], barrier_array_2, dot_ports_2xN[1, i], _ = dot_creation(gate_length_addition[i], 
                            dot_positions_2xN[1, i], dot_radius, angle_resolution,
                            layer_dot[layer], dot_width_gate, dot_tilts_left_2xN[0, i],
                            dot_tilts_right_2xN[0, i], dot_rotations_2xN[1, i], layer_barrier[layer], [], [])
        D << dot_2xN[0, i]
        D << dot_gates[0, i]
        D << barrier_array_1
        D << dot_2xN[1, i]
        D << dot_gates[1, i]
        D << barrier_array_2
        
    """
    Build the sensor dots, and their barriers, based on the positions and parameters defined above.
    """
    
    sensor_dot_gate_length = [300, 300]
    ports = [['1', '3'], ['3', '1']]
    port_routings = [[[50, 70], [70, 50]], [[70, 50], [50, 70]]]
    
    for i in range(len(sensordot_positions_2xN)):
        sensor_dots_2xN[i], gate, barrier_array, sensor_dot_ports_2xN[i], touple_barrier_ports = dot_creation(sensor_dot_gate_length[i],
                            sensordot_positions_2xN[i], sensordot_radius, angle_resolution,
                            layer_dot[layer], sensor_width_gate, sensor_tilts_left_2x2[i],
                            sensor_tilts_right_2x2[i], sensor_rotations_2x2[i], layer_barrier[layer],
                            ports[i], port_routings[i], barrier_rotations_2x2, sensor_routing_directions[i])
        D << sensor_dots_2xN[i]
        D << gate
        D << barrier_array
        embedded_barrier_ports_2x2.append(touple_barrier_ports[0])
        embedded_barrier_ports_2x2.append(touple_barrier_ports[1])
    
    
    """
    Manually create the screening gates, and their ports. The number of screening gates
    is fixed, but some of the parameters are not, such as the length of it.
    """
    
    center_dots_screening_1 = pg.polygon_ports(xpts =[dot_gates[1, 0].xmin - overhang_screening,
                                                      dot_gates[1, amount_of_dots - 1].xmax + overhang_screening, 
                                                      dot_gates[1, amount_of_dots - 1].xmax + overhang_screening, 
                                                      dot_gates[1, 0].xmin - overhang_screening], 
                               ypts = [               dot_2xN[1, 0].center[1] + dot_radius + length_dot_perimeter_to_screening, 
                                                      dot_2xN[1, 0].center[1] + dot_radius + length_dot_perimeter_to_screening, 
                                                      dot_2xN[1, 0].center[1] + dot_radius + screen_length + length_dot_perimeter_to_screening,
                                                      dot_2xN[1, 0].center[1] + dot_radius + screen_length + length_dot_perimeter_to_screening],
                               layer = layer_screening[layer])
    
    screening_gate_1 = pg.polygon_ports(xpts =  [dot_2xN[1, 0].center[0] + dot_radius - 0.65*screening_gate_width,
                                                 dot_2xN[1, 0].center[0] + dot_radius + 0.3*screening_gate_width, 
                                                 dot_2xN[1, 0].center[0] + dot_radius + 0.3*screening_gate_width, 
                                                 dot_2xN[1, 0].center[0] + dot_radius - 0.65*screening_gate_width], 
                               ypts = [          0.95*center_dots_screening_1.ymax, 
                                                 0.95*center_dots_screening_1.ymax, 
                                                 0.95*center_dots_screening_1.ymax + gate_length_screen + 30,
                                                 0.95*center_dots_screening_1.ymax + gate_length_screen],
                               layer = layer_screening[layer])


    center_dots_screening_2 = pg.polygon_ports(xpts =[dot_gates[0, 0].xmin - overhang_screening,
                                                      dot_gates[0, amount_of_dots - 1].xmax + overhang_screening, 
                                                      dot_gates[0, amount_of_dots - 1].xmax + overhang_screening, 
                                                      dot_gates[0, 0].xmin - overhang_screening], 
                               ypts = [               dot_2xN[0, 0].center[1] - dot_radius - length_dot_perimeter_to_screening, 
                                                      dot_2xN[0, 0].center[1] - dot_radius - length_dot_perimeter_to_screening, 
                                                      dot_2xN[0, 0].center[1] - dot_radius - screen_length - length_dot_perimeter_to_screening,
                                                      dot_2xN[0, 0].center[1] - dot_radius - screen_length - length_dot_perimeter_to_screening],
                               layer = layer_screening[layer])
    
    screening_gate_2 = pg.polygon_ports(xpts =  [dot_2xN[0, amount_of_dots - 1].center[0] - dot_radius - 0.3*screening_gate_width,
                                                 dot_2xN[0, amount_of_dots - 1].center[0] - dot_radius + 0.65*screening_gate_width, 
                                                 dot_2xN[0, amount_of_dots - 1].center[0] - dot_radius + 0.65*screening_gate_width, 
                                                 dot_2xN[0, amount_of_dots - 1].center[0] - dot_radius - 0.3*screening_gate_width], 
                               ypts = [          0.95*center_dots_screening_2.ymin, 
                                                 0.95*center_dots_screening_2.ymin, 
                                                 0.95*center_dots_screening_2.ymin - gate_length_screen,
                                                 0.95*center_dots_screening_2.ymin - gate_length_screen - 30],
                               layer = layer_screening[layer])
    
    sensor_dot_screening_1 = pg.polygon_ports(xpts =    [sensor_dots_2xN[0].xmin - distance_from_sensordot_to_screening_gate - 50,
                                                         sensor_dots_2xN[0].xmin - distance_from_sensordot_to_screening_gate, 
                                                         sensor_dots_2xN[0].xmin - distance_from_sensordot_to_screening_gate, 
                                                         sensor_dots_2xN[0].xmin - distance_from_sensordot_to_screening_gate - 70], 
                               ypts = [                  sensor_dots_2xN[0].center[1] - dot_radius, 
                                                         sensor_dots_2xN[0].center[1] - dot_radius, 
                                                         sensor_dots_2xN[0].center[1] + 1.3*dot_radius,
                                                         sensor_dots_2xN[0].center[1] + 0.7*dot_radius], layer = layer_screening[layer])
    
    sensor_dot_screening_2 = pg.polygon_ports(xpts =    [sensor_dots_2xN[1].xmax + distance_from_sensordot_to_screening_gate + 50,
                                                         sensor_dots_2xN[1].xmax + distance_from_sensordot_to_screening_gate, 
                                                         sensor_dots_2xN[1].xmax + distance_from_sensordot_to_screening_gate, 
                                                         sensor_dots_2xN[1].xmax + distance_from_sensordot_to_screening_gate + 70], 
                               ypts = [                  sensor_dots_2xN[1].center[1] - dot_radius, 
                                                         sensor_dots_2xN[1].center[1] - dot_radius, 
                                                         sensor_dots_2xN[1].center[1] + 1.3*dot_radius,
                                                         sensor_dots_2xN[1].center[1] + 0.7*dot_radius], layer = layer_screening[layer])
    
    screening_gates = [center_dots_screening_1, screening_gate_1, center_dots_screening_2, screening_gate_2,
                       sensor_dot_screening_1, sensor_dot_screening_2]
    
    for i in screening_gates:
        D << i                  # Add the screening gates to the device
    
    screening_ports = []

    # Add the screening gate ports to the screening_ports array, such that they can be routed out.
    
    screening_ports.append(sensor_dot_screening_1.ports['3'])
    screening_ports.append(sensor_dot_screening_2.ports['3'])
    screening_ports.append(screening_gate_1.ports['3'])
    screening_ports.append(screening_gate_2.ports['3'])
    
    
    """
    Creating the arrays and such for the inner barrier crosses. 
    """
    
    
    upper_lower_cross_positions = np.zeros((2, amount_of_dots - 1), dtype = tuple)
    middle_cross_positions = np.zeros(amount_of_dots, dtype = tuple)
        
    upper_lower_crosses = np.zeros((2, amount_of_dots - 1), dtype = tuple)
    middle_crosses = np.zeros(amount_of_dots, dtype = tuple)
    
    upper_lower_cross_ports = np.zeros((2, amount_of_dots - 1), dtype = tuple)
    middle_cross_ports = np.zeros(amount_of_dots, dtype = tuple)
    
    upper_lower_cross_direction = np.zeros((2, amount_of_dots - 1), dtype = str)
    middle_cross_directions = np.zeros(amount_of_dots, dtype = str)
    
    upper_lower_routing_direction = np.zeros((2, amount_of_dots - 1), dtype = tuple)
    middle_routing_direction = np.zeros(amount_of_dots, dtype = tuple)
    
    routing_middle_temp = np.asarray(['1', '3']*amount_of_dots) # Make an alternating array of 1's and 3's, so that the routing for the inner crosses alternates between up and down
    
    middle_cross_directions[0] = '3' # Manually set the first two directions, since they anchored due to the sensor dots
    
    middle_cross_directions[amount_of_dots - 1] = '1'
    
    """ 
    The code directly below sets all the middle routing directions for the 
    cross barriers in the center.
    Below that, we set the positions of the upper and lower cross barriers, and we
    start setting their routing direction. On the end, middle cross barriers,
    we also set a skew on their ports, so no thinning of the gates occur.
    """
    
    for i in range(1, amount_of_dots - 1):
        middle_cross_directions[i] = routing_middle_temp[i]
    
    for i in range(amount_of_dots):
        if i < amount_of_dots - 1:
            upper_lower_cross_positions[0, i] = [dot_2xN[0, i].center[0] + 0.5*spacing_dots,
                                                 dot_2xN[0, i].center[1]]
            upper_lower_cross_positions[1, i] = [dot_2xN[1, i].center[0] + 0.5*spacing_dots,
                                                 dot_2xN[1, i].center[1]]
            middle_cross_positions[i] = [dot_2xN[0, i].center[0], dot_2xN[0, i].center[1] + 0.5*spacing_dots]
            upper_lower_cross_direction[0, i] = '1'
            upper_lower_cross_direction[1, i] = '3'
            upper_lower_routing_direction[0, i] = [0, 0]
            upper_lower_routing_direction[1, i] = [0, 0]
            middle_routing_direction[i] = [0, 0]
        if i == 0:
            middle_routing_direction[i] = [0, 20]
        if i == amount_of_dots - 1:
            middle_cross_positions[i] = [dot_2xN[0, i].center[0], dot_2xN[0, i].center[1] + 0.5*spacing_dots]
            middle_routing_direction[i] = [0, 20]

    
    """ 
    Going through the indeces 0 and 1, to indicate lower and upper crosses.
    The middle cross barriers are only created when the loop is at index 0, since
    there is only one row of them. 
    """
    
    for j in range(2):
        for i in range(amount_of_dots):
            if i < amount_of_dots - 1:
                cross_barrier_obj, cross_barrier_port = cross_barriers(mid_point = upper_lower_cross_positions[j, i],
                                                                    width = width_cross,
                                                                    length = length_cross, layer = layer_barrier[layer],
                                                                    port_cardinal = upper_lower_cross_direction[j, i],
                                                                    port_routing_direction = upper_lower_routing_direction[j, i])
                upper_lower_crosses[j, i] = cross_barrier_obj
                appendage, appendage_port = make_appendage(cross_barrier_obj, upper_lower_cross_direction[j, i],
                                                           appendage_length_addition[i], layer_barrier, layer)
                upper_lower_cross_ports[j, i] = appendage_port
                if j == 0:
                    cross_barrier_obj_middle, cross_barrier_port_middle = cross_barriers(
                                                                        mid_point = middle_cross_positions[i],
                                                                        width = width_cross,
                                                                        length = length_cross, layer = layer_barrier[layer],
                                                                        port_cardinal = middle_cross_directions[i],
                                                                        port_routing_direction = middle_routing_direction[i])
                    middle_crosses[i] = cross_barrier_obj_middle
                    if i != 0 and i != amount_of_dots - 1:
                        appendage_middle, appendage_port_middle = make_appendage(cross_barrier_obj_middle, middle_cross_directions[i], 
                                                                                 appendage_length_middle_addition[i], layer_barrier, layer)
                        middle_cross_ports[i] = appendage_port_middle
                        D << appendage_middle
                    else: 
                        middle_cross_ports[i] = cross_barrier_port_middle
                D << middle_crosses[i]
                D << upper_lower_crosses[j, i]
                D << appendage
            if i == amount_of_dots - 1 and j == 0:
                cross_barrier_obj_middle, cross_barrier_port_middle = cross_barriers(
                                                                        mid_point = middle_cross_positions[i],
                                                                        width = width_cross,
                                                                        length = length_cross, layer = layer_barrier[layer],
                                                                        port_cardinal = middle_cross_directions[i],
                                                                        port_routing_direction = middle_routing_direction[i])
                middle_crosses[i] = cross_barrier_obj_middle
                middle_cross_ports[i] = cross_barrier_port_middle
                D << middle_crosses[i]
    
    """
    Build the barriers in between the sensors dots and the dots and add their ports
    to a list.
    """
    
    barriers_between_sensors_and_dots = np.zeros(2, dtype = tuple)
    
    # barriers_between_sensors_and_dots[0] = pg.polygon_ports(xpts = [(dot_2xN[0, 0].center[0] - abs(sensor_dots_2xN[0].center[0]))/2 - 0.5*width_barrier + 0.5*sensor_dot_dot_radii_difference,
    #                                                    (dot_2xN[0, 0].center[0] - abs(sensor_dots_2xN[0].center[0]))/2 + 0.5*width_barrier + 0.5*sensor_dot_dot_radii_difference, 
    #                                                    (dot_2xN[0, 0].center[0] - abs(sensor_dots_2xN[0].center[0]))/2 + 0.5*width_barrier + 0.5*sensor_dot_dot_radii_difference, 
    #                                                    (dot_2xN[0, 0].center[0] - abs(sensor_dots_2xN[0].center[0]))/2 - 0.5*width_barrier + 0.5*sensor_dot_dot_radii_difference], 
    #                            ypts = [                sensor_dots_2xN[0].center[1] - lenght_barrier_between_sensor_and_dot, 
    #                                                    sensor_dots_2xN[0].center[1] - lenght_barrier_between_sensor_and_dot - tilt_barrier_between_sensor_and_dot, 
    #                                                    sensor_dots_2xN[0].center[1] + 0.7*dot_radius,
    #                                                    sensor_dots_2xN[0].center[1] + 0.7*dot_radius], layer = layer_barrier[layer])
    
    vert_hori_1, barriers_between_sensors_and_dots[0] = cross_barriers(mid_point = [(dot_2xN[0, 0].center[0] - abs(sensor_dots_2xN[0].center[0]))/2 + (sensor_dot_dot_radii_difference)/2, sensor_dots_2xN[0].center[1]],
                                                                    width = width_cross,
                                                                    length = length_cross, layer = layer_barrier[layer],
                                                                    port_cardinal = '1',
                                                                    port_routing_direction = [200, 180])
    
    
    # barriers_between_sensors_and_dots[1] = pg.polygon_ports(xpts = [(dot_2xN[1, amount_of_dots - 1].center[0] + abs(sensor_dots_2xN[1].center[0]))/2 - 0.5*width_barrier - 0.5*sensor_dot_dot_radii_difference,
    #                                                    (dot_2xN[1, amount_of_dots - 1].center[0] + abs(sensor_dots_2xN[1].center[0]))/2 + 0.5*width_barrier - 0.5*sensor_dot_dot_radii_difference, 
    #                                                    (dot_2xN[1, amount_of_dots - 1].center[0] + abs(sensor_dots_2xN[1].center[0]))/2 + 0.5*width_barrier - 0.5*sensor_dot_dot_radii_difference, 
    #                                                    (dot_2xN[1, amount_of_dots - 1].center[0] + abs(sensor_dots_2xN[1].center[0]))/2 - 0.5*width_barrier - 0.5*sensor_dot_dot_radii_difference], 
    #                            ypts = [                sensor_dots_2xN[1].center[1] - 0.7*dot_radius, 
    #                                                    sensor_dots_2xN[1].center[1] - 0.7*dot_radius, 
    #                                                    sensor_dots_2xN[1].center[1] + lenght_barrier_between_sensor_and_dot,
    #                                                    sensor_dots_2xN[1].center[1] + lenght_barrier_between_sensor_and_dot + tilt_barrier_between_sensor_and_dot], layer = layer_barrier[layer])
    
    vert_hori_2, barriers_between_sensors_and_dots[1] = cross_barriers(mid_point = [sensor_dots_2xN[1].center[0] - (sensor_dots_2xN[1].center[0]- dot_2xN[1, -1].center[0])/2 - sensor_dot_dot_radii_difference/2, sensor_dots_2xN[1].center[1]],
                                                                width = width_cross,
                                                                length = length_cross, layer = layer_barrier[layer],
                                                                port_cardinal = '3',
                                                                port_routing_direction = [200, 180])
    
    barrier_ports = []
    
    barrier_ports.append(vert_hori_1[1].ports['1'])
    barrier_ports.append(vert_hori_2[1].ports['3'])
    
    D << vert_hori_1[0]
    D << vert_hori_2[0]
    D << vert_hori_1[1]
    D << vert_hori_2[1]
    
    for i in embedded_barrier_ports_2x2:
        barrier_ports.append(i)
    
    # for objec in barriers_between_sensors_and_dots:
    #     D << objec
    
    """ 
    Build the ohmics and their ports and add them to a list. Again,
    in this design, the ohmics do not vary with the amount of dots, 
    so they can be created manually.
    """
    
    distance_from_centersensordot_to_ohmic = sensordot_radius + dot_params['distance_from_sensordot_to_ohmic']
    
    distance_from_sensordot_to_ohmic = dot_params['distance_from_sensordot_to_ohmic']
    
    ohmics = np.zeros(4, dtype = tuple)
    
    ohmics[0] = pg.polygon_ports(xpts = [            sensor_dots_2xN[0].center[0] - 0.35*sensordot_radius,
                                                     sensor_dots_2xN[0].center[0] + 0.35*sensordot_radius , 
                                                     sensor_dots_2xN[0].center[0] + 0.35*sensordot_radius, 
                                                     sensor_dots_2xN[0].center[0] - 0.4*sensordot_radius], 
                               ypts = [              sensor_dots_2xN[0].center[1] - distance_from_centersensordot_to_ohmic, 
                                                     sensor_dots_2xN[0].center[1] - distance_from_centersensordot_to_ohmic, 
                                                     sensor_dots_2xN[0].center[1] - 2.3*sensordot_radius - distance_from_sensordot_to_ohmic,
                                                     sensor_dots_2xN[0].center[1] - 1.8*sensordot_radius - distance_from_sensordot_to_ohmic],
                               layer = layer_ohmic[layer])
    
    ohmics[1] = pg.polygon_ports(xpts = [            sensor_dots_2xN[0].center[0] - 0.4*sensordot_radius,
                                                     sensor_dots_2xN[0].center[0] + 0.35*sensordot_radius , 
                                                     sensor_dots_2xN[0].center[0] + 0.35*sensordot_radius, 
                                                     sensor_dots_2xN[0].center[0] - 0.35*sensordot_radius], 
                               ypts = [              sensor_dots_2xN[0].center[1] + 1.8*sensordot_radius + distance_from_sensordot_to_ohmic, 
                                                     sensor_dots_2xN[0].center[1] + 2.3*sensordot_radius + distance_from_sensordot_to_ohmic, 
                                                     sensor_dots_2xN[0].center[1] + distance_from_centersensordot_to_ohmic,
                                                     sensor_dots_2xN[0].center[1] + distance_from_centersensordot_to_ohmic],
                               layer = layer_ohmic[layer])
    
    ohmics[2] = pg.polygon_ports(xpts = [            sensor_dots_2xN[1].center[0] - 0.35*sensordot_radius,
                                                     sensor_dots_2xN[1].center[0] - 0.35*sensordot_radius , 
                                                     sensor_dots_2xN[1].center[0] + 0.4*sensordot_radius, 
                                                     sensor_dots_2xN[1].center[0] + 0.35*sensordot_radius], 
                               ypts = [              sensor_dots_2xN[1].center[1] - distance_from_centersensordot_to_ohmic, 
                                                     sensor_dots_2xN[1].center[1] - 2.3*sensordot_radius - distance_from_sensordot_to_ohmic, 
                                                     sensor_dots_2xN[1].center[1] - 1.8*sensordot_radius - distance_from_sensordot_to_ohmic,
                                                     sensor_dots_2xN[1].center[1] - distance_from_centersensordot_to_ohmic],
                               layer = layer_ohmic[layer])
    
    ohmics[3] = pg.polygon_ports(xpts = [            sensor_dots_2xN[1].center[0] - 0.35*sensordot_radius,
                                                     sensor_dots_2xN[1].center[0] + 0.35*sensordot_radius , 
                                                     sensor_dots_2xN[1].center[0] + 0.4*sensordot_radius, 
                                                     sensor_dots_2xN[1].center[0] - 0.35*sensordot_radius], 
                               ypts = [              sensor_dots_2xN[1].center[1] + distance_from_centersensordot_to_ohmic, 
                                                     sensor_dots_2xN[1].center[1] + distance_from_centersensordot_to_ohmic, 
                                                     sensor_dots_2xN[1].center[1] + 1.8*sensordot_radius + distance_from_sensordot_to_ohmic,
                                                     sensor_dots_2xN[1].center[1] + 2.3*sensordot_radius + distance_from_sensordot_to_ohmic],
                               layer = layer_ohmic[layer])
    
    for objec in ohmics:
        D << objec
        
    ohmic_ports = [ohmics[0].ports['3'], ohmics[1].ports['1'], ohmics[2].ports['2'], ohmics[3].ports['3']]
    
    ext_ports = [[], [], [], []]
    
    """ 
    Setting the actual placements for the ports.
    """
    
    North_placements, South_placements, West_placements, East_placements = calculate_placements(amount_of_dots)
    N_bond_numbers, S_bond_numbers, W_bond_numbers, E_bond_numbers = calculate_placements(amount_of_dots)
        
    North_length = len(North_placements)            # Having calculated the number of routings, 
                                                    # we use this number for all
    South_length = len(South_placements)            # subsequent calculations.

    ext_port_N, ext_port_S, ext_port_W, ext_port_E = [], [], [], []
    
    bond_numbers = [N_bond_numbers, S_bond_numbers, W_bond_numbers, E_bond_numbers]

    list_of_numbers = [North_placements, South_placements, West_placements, East_placements]

    """ 
    Going through each list of numbers and adding the external ports.
    """

    counter_1 = 0

    distance_ext_port_north = system_params['distance_external_port_north']*1000
    distance_ext_port_south = system_params['distance_external_port_south']*1000
    distance_ext_port_west = system_params['distance_external_port_west']*1000
    distance_ext_port_east = system_params['distance_external_port_east']*1000

    for i in list_of_numbers:
        for j in range(len(i)):
            if counter_1 == 0:
                ext_port_N.append('EXT_N' + str(j))
                ext_port_N[j] = D.add_port(name='EXT_N' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]), distance_ext_port_north], width = system_params['width_external_port'], orientation = 270)
            if counter_1 == 1:
                ext_port_S.append('EXT_S' + str(j))
                ext_port_S[j] = D.add_port(name='EXT_S' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]), distance_ext_port_south], width = system_params['width_external_port'], orientation = 90)
            if counter_1 == 2:
                ext_port_W.append('EXT_W' + str(j))
                ext_port_W[j] = D.add_port(name='EXT_W' + str(j), midpoint = [distance_ext_port_west, y_midpoint_ext_port_v(i[j])], width = system_params['width_external_port'], orientation = 0)
            if counter_1 == 3:
                ext_port_E.append('EXT_E' + str(j))
                ext_port_E[j] = D.add_port(name='EXT_E' + str(j), midpoint = [distance_ext_port_east, y_midpoint_ext_port_v(i[j])], width = system_params['width_external_port'], orientation = 180)
        counter_1 += 1
    
    """ 
    Each routed port and adjoining gate needs a layer, so those
    are created below.
    """
    
    layer_N = np.zeros(North_length, dtype = tuple)
    layer_S = np.zeros(South_length, dtype = tuple)
    
    layer_N_pads = np.zeros(North_length, dtype = tuple)
    layer_S_pads = np.zeros(South_length, dtype = tuple)
    layer_W_pads = np.zeros(North_length, dtype = tuple)
    layer_E_pads = np.zeros(South_length, dtype = tuple)
    
    layer_N[0] = layer_barrier[layer]
    layer_N[1] = layer_dot[layer]
    layer_N[2] = layer_screening[layer]
    
    """
    The current set-up for the outer layers is as follows:
        layers_for_outer[0] = dot
        layers_for_outer[1] = barrier
        layers_for_outer[2] = screening
        layers_for_outer[3] = ohmic 
    The outer layers is meant for the creation of the outermost gates,
    in each of the four designs that fit on a chip.
    This is done so that none of them overlap with layers
    from the inner gates.
    """
    
    layer_N_pads[0] = layers_for_outer[1]
    layer_N_pads[1] = layers_for_outer[0]
    layer_N_pads[2] = layers_for_outer[2]
    
    layer_N[North_length - 1] = layer_barrier[layer]
    
    layer_N_pads[North_length - 1] = layers_for_outer[1]
    
    
    layer_S[0] = layer_barrier[layer]
    layer_S[1] = layer_dot[layer]
    
    layer_S_pads[0] = layers_for_outer[1]
    layer_S_pads[1] = layers_for_outer[0]
    
    """ 
    Initialising the list of layers that we will draw from. This is the general 
    form of the layers, so we simply multiply it by some integer to get a long, 
    repeating list.
    """
    
    North_layers_list = [layer_barrier[layer], layer_dot[layer], layer_barrier[layer], layer_barrier[layer], layer_dot[layer]]*20
    
    North_layers_pad_list = [layers_for_outer[1], layers_for_outer[0], layers_for_outer[1], layers_for_outer[1], layers_for_outer[0]]*20
    
    South_layers_list = [layer_barrier[layer], layer_dot[layer], layer_barrier[layer], layer_barrier[layer], layer_dot[layer]]*20

    South_layers_pad_list = [layers_for_outer[1], layers_for_outer[0], layers_for_outer[1], layers_for_outer[1], layers_for_outer[0]]*20
    
    """ 
    Adding the various layers to each array in the correct order. 
    Start by initialising each end of the arrays,
    as those to not vary based on the amount of dots.
    
    """
    
    if amount_of_dots > 1:
        layer_N[North_length - 2] = layer_dot[layer]
        layer_S[South_length - 1] = layer_barrier[layer]
        layer_S[South_length - 2] = layer_dot[layer]
        layer_S[South_length - 3] = layer_screening[layer]
        
        layer_N_pads[North_length - 2] = layers_for_outer[0]
        layer_S_pads[South_length - 1] = layers_for_outer[1]
        layer_S_pads[South_length - 2] = layers_for_outer[0]
        layer_S_pads[South_length - 3] = layers_for_outer[2]
        counter = 0
        for i in range(3, North_length - 2):
            layer_N[i] = North_layers_list[counter]
            layer_N_pads[i] = North_layers_pad_list[counter]
            counter += 1
        counter = 0
        for i in range(2, South_length - 3):
            layer_S[i] = South_layers_list[counter]
            layer_S_pads[i] = South_layers_pad_list[counter]
            counter += 1
    else:
        layer_S[South_length - 1] = layer_screening[layer]
        layer_S_pads[South_length - 1] = layers_for_outer[2]
    
    layer_W = [layer_ohmic[layer], layer_barrier[layer], layer_dot[layer], layer_screening[layer], layer_barrier[layer], layer_ohmic[layer]]
    layer_E = [layer_ohmic[layer], layer_barrier[layer], layer_dot[layer], layer_screening[layer], layer_barrier[layer], layer_ohmic[layer]]

    layer_W_pads = [layers_for_outer[3], layers_for_outer[1], layers_for_outer[0], layers_for_outer[2], layers_for_outer[1], layers_for_outer[3]]
    layer_E_pads = [layers_for_outer[3], layers_for_outer[1], layers_for_outer[0], layers_for_outer[2], layers_for_outer[1], layers_for_outer[3]]

    if amount_of_dots == 1:
        North_objects_to_route = np.zeros(North_length - 1, dtype = tuple)
    
        South_objects_to_route = np.zeros(South_length - 1, dtype = tuple)
    else:
        North_objects_to_route = np.zeros(North_length, dtype = tuple)
    
        South_objects_to_route = np.zeros(South_length, dtype = tuple)
    
    North_objects_to_route[0], North_objects_to_route[1], North_objects_to_route[2] = middle_cross_ports[0], dot_ports_2xN[1, 0], screening_ports[2]
    
    South_objects_to_route[0], South_objects_to_route[1] = barrier_ports[0], dot_ports_2xN[0, 0]
    
    North_objects_to_route[len(North_objects_to_route) - 1] = barrier_ports[1]

    West_objects_to_route = [ohmic_ports[0], barrier_ports[2], sensor_dot_ports_2xN[0], screening_ports[0],
                              barrier_ports[3], ohmic_ports[1]]
    
    East_objects_to_route = [ohmic_ports[2], barrier_ports[4], sensor_dot_ports_2xN[1], screening_ports[1],
                             barrier_ports[5], ohmic_ports[3]]


    def port_lists(multiple):
        
        """
        This function creates lists of ports for the North and South sides of the device.
        It works by looking at mupltiples of 5, since that is the repeating
        pattern I found in the current design.
        As it is, it works for a smaller amount of dots, but more functionalty is needed
        for amount_of_dots > 9.
        """
        
        
        end = multiple*5
        North_ports_list = np.zeros(North_length, dtype = tuple)
        South_ports_list = np.zeros(South_length, dtype = tuple)
        upper_lower_1_counter, upper_lower_0_counter, middle_cross_counter, dot_ports_0_2xN_counter = 0, 0, 1, 1
        dot_ports_1_2xN_counter = 1
        for i in range(0, end, 5):
            if upper_lower_1_counter < len(upper_lower_cross_ports[1]):
                North_ports_list[i] = upper_lower_cross_ports[1, upper_lower_1_counter]
                upper_lower_1_counter += 1
                
            if upper_lower_0_counter < len(upper_lower_cross_ports[0]):
                South_ports_list[i] = upper_lower_cross_ports[0, upper_lower_0_counter]
                upper_lower_0_counter += 1
                
            if dot_ports_1_2xN_counter < len(dot_ports_2xN[1]):
                North_ports_list[i + 1] = dot_ports_2xN[1, dot_ports_1_2xN_counter]
                dot_ports_1_2xN_counter += 1
                
            if dot_ports_0_2xN_counter < len(dot_ports_2xN[1]):
                South_ports_list[i + 1] = dot_ports_2xN[0, dot_ports_0_2xN_counter]
                dot_ports_0_2xN_counter += 1
            
            if middle_cross_counter < len(middle_cross_ports):
                North_ports_list[i + 2] = middle_cross_ports[middle_cross_counter]
                middle_cross_counter += 1
                
            if upper_lower_0_counter < len(upper_lower_cross_ports[0]):
                South_ports_list[i + 2] = upper_lower_cross_ports[0, upper_lower_0_counter]
                upper_lower_0_counter += 1
                
            if upper_lower_1_counter < len(upper_lower_cross_ports[1]):
                North_ports_list[i + 3] = upper_lower_cross_ports[1, upper_lower_1_counter]
                upper_lower_1_counter += 1
                
            if middle_cross_counter < len(middle_cross_ports):
                South_ports_list[i + 3] = middle_cross_ports[middle_cross_counter]
                middle_cross_counter += 1
                
            if dot_ports_1_2xN_counter < len(dot_ports_2xN[1]):
                North_ports_list[i + 4] = dot_ports_2xN[1, dot_ports_1_2xN_counter]
                dot_ports_1_2xN_counter += 1
                
            if dot_ports_0_2xN_counter < len(dot_ports_2xN[0]):
                South_ports_list[i + 4] = dot_ports_2xN[0, dot_ports_0_2xN_counter]
                dot_ports_0_2xN_counter += 1
                

        return North_ports_list, South_ports_list

    """ 
    Creating the porting list and filling it. 
    """

    if amount_of_dots >= 2:
        multiple = 2

    if amount_of_dots >= 4:
        multiple = 4
    
    if amount_of_dots >= 6:
        multiple = 6
    
    if amount_of_dots >= 8:
        multiple = 8
    
    if amount_of_dots >= 10:
        multiple = 10
    

    North_ports_list, South_ports_list = port_lists(multiple)
    
    """ 
    Setting the end objects, which do not vary with the amount of dots.
    These are the 2-3 first objects on each end.
    """
    
    if amount_of_dots > 1:
        North_objects_to_route[North_length - 2] = dot_ports_2xN[1, amount_of_dots - 1]
        South_objects_to_route[South_length - 1] = middle_cross_ports[len(middle_cross_ports) - 1]
        South_objects_to_route[South_length - 2] = dot_ports_2xN[0, amount_of_dots - 1]
        South_objects_to_route[South_length - 3] = screening_ports[3]
        counter = 0
        for i in range(3, North_length - 2):
            North_objects_to_route[i] = North_ports_list[counter]
            counter += 1
        counter = 0
        for i in range(2, South_length - 3):
            South_objects_to_route[i] = South_ports_list[counter]
            counter += 1
    else:
        South_objects_to_route[len(South_objects_to_route) - 1] = screening_ports[3]
    
    """ 
    Creating the arrays of external ports to route to, filling it up based
    on the above defined objects, and finally routing them.
    """

    if amount_of_dots == 1:
        ext_ports = [ext_port_N[:len(North_objects_to_route) - 1],
                     ext_port_S[:len(South_objects_to_route) - 1], ext_port_W, ext_port_E]
    else:
        ext_ports = [ext_port_N, ext_port_S, ext_port_W, ext_port_E]

    route_ext_N = []
    route_ext_S = []
    route_ext_W = []
    route_ext_E = []


    for i in ext_ports:
        for j in range(len(i)):
            if i == ext_port_N:
                route_ext_N.append('route_ext_N' + str(j))
                route_ext_N[j] = D.add_ref(pr.route_quad(North_objects_to_route[j], i[j], layer = layer_N[j]))
            if i == ext_port_S:
                route_ext_S.append('route_ext_S' + str(j))
                route_ext_S[j] = D.add_ref(pr.route_quad(South_objects_to_route[j], i[j], layer = layer_S[j]))
            if i == ext_port_W:
                route_ext_W.append('route_ext_W' + str(j))
                route_ext_W[j] = D.add_ref(pr.route_quad(West_objects_to_route[j], i[j], layer = layer_W[j]))
            if i == ext_port_E:
                route_ext_E.append('route_ext_E' + str(j))
                route_ext_E[j] = D.add_ref(pr.route_quad(East_objects_to_route[j], i[j], layer = layer_E[j]))

    """ 
    Making a new device, without markers, so that it can be output at the end of the
    function.
    """

    D_without_local_markers = pg.union(D, by_layer = True)

    D_local_markers = Device()

    d_ref6 = D_local_markers.add_ref(D)

    # bond_extentable = True                  # Boolean which decides whether or not
    #                                         # we constrain the x-length

    
    bond_extendable = system_params['bond_extendable']

    """ 
    The code below is for added funtionality, in which the placement of the inner
    markers are dependent on whether or not we have constrained the x-length of
    each design. This is not generally meant for the completed designs, but for 
    creating new designs.
    """
    
    D_local_markers = Device()              # Making the structure with local markers

    D_local_markers.add_ref(D)
    
    nr_of_locals = len(local_markers_pos)
    
    inner_NE = []
    inner_NW = []
    inner_SE = []
    inner_SW = []
    
    local_marker_positions_NE = []
    local_marker_positions_NW = []
    local_marker_positions_SE = []
    local_marker_positions_SW = []
    
    multipliers = [40000, 50000, 60000, 70000]
    
    for i in range(nr_of_locals):
        local_marker_positions_NE.append([D_local_markers.center[0] + math.sqrt(1/2)*local_markers_pos[i]*1000, D_local_markers.center[1] + math.sqrt(1/2)*local_markers_pos[i]*1000 - multipliers[i]])
        local_marker_positions_NW.append([D_local_markers.center[0] - math.sqrt(1/2)*local_markers_pos[i]*1000, D_local_markers.center[1] + math.sqrt(1/2)*local_markers_pos[i]*1000 - multipliers[i]])
        local_marker_positions_SE.append([D_local_markers.center[0] + math.sqrt(1/2)*local_markers_pos[i]*1000, D_local_markers.center[1] - math.sqrt(1/2)*local_markers_pos[i]*1000 + multipliers[i]])
        local_marker_positions_SW.append([D_local_markers.center[0] - math.sqrt(1/2)*local_markers_pos[i]*1000, D_local_markers.center[1] - math.sqrt(1/2)*local_markers_pos[i]*1000 + multipliers[i]])

    if bond_extendable == True:
        local_marker_positions = {}
        for i in range(nr_of_locals):
            inner_SW.append("inner_marker_SW" + str(i))
            inner_SW[i] = inner_marker_generator(marker_params['width_local_marker']*1000, marker_params['length_local_marker']*1000, layer_marker, cardinal = "A" + str(i))
            inner_SW[i].move(destination = (local_marker_positions_SW[i]))
            D_local_markers.add_ref(inner_SW[i], alias = "inner_marker_SW" + str(i))
            local_marker_positions.update({"A" + str(i): [inner_SW[i].center[0], inner_SW[0].center[1]]})
            
            inner_SE.append("inner_marker_SE" + str(i))
            inner_SE[i] = inner_marker_generator(marker_params['width_local_marker']*1000, marker_params['length_local_marker']*1000, layer_marker, cardinal = "B" + str(i))
            inner_SE[i].move(destination = (local_marker_positions_SE[i]))
            D_local_markers.add_ref(inner_SE[i], alias = "inner_marker_SE" + str(i))
            local_marker_positions.update({"B" + str(i): [inner_SE[i].center[0], inner_SE[0].center[1]]})
            
            inner_NW.append("inner_marker_NW" + str(i))
            inner_NW[i] = inner_marker_generator(marker_params['width_local_marker']*1000, marker_params['length_local_marker']*1000, layer_marker, cardinal = "C" + str(i))
            inner_NW[i].move(destination = (local_marker_positions_NW[i]))
            D_local_markers.add_ref(inner_NW[i], alias = "inner_marker_NW" + str(i))
            local_marker_positions.update({"C" + str(i): [inner_NW[i].center[0], inner_NW[0].center[1]]})
            
            inner_NE.append("inner_marker_NE" + str(i))
            inner_NE[i] = inner_marker_generator(marker_params['width_local_marker']*1000, marker_params['length_local_marker']*1000, layer_marker, cardinal = "D" + str(i))
            inner_NE[i].move(destination = (local_marker_positions_NE[i]))
            D_local_markers.add_ref(inner_NE[i], alias = "inner_marker_NE" + str(i))
            local_marker_positions.update({"D" + str(i): [inner_NE[i].center[0], inner_NE[0].center[1]]})
    else:
        inner_marker_NE = inner_marker_generator(system_params['width_local_marker']*1000, system_params['length_local_marker']*1000, layer_marker, cardinal = 'D')
        inner_marker_NE.move(destination = (np.cos(45)*system_params['length_to_local']*1000, i*0.825*np.cos(45)*system_params['length_to_local']*1000))
        D_local_markers.add_ref(inner_marker_NE, alias = "inner_marker_NE")

        inner_marker_NW = inner_marker_generator(system_params['width_local_marker']*1000, system_params['length_local_marker']*1000, layer_marker, cardinal = 'C')
        inner_marker_NW.move(destination = (-np.cos(45)*system_params['length_to_local']*1000, 0.825*np.cos(45)*system_params['length_to_local']*1000))
        D_local_markers.add_ref(inner_marker_SE, alias = "inner_marker_SE")

        inner_marker_SE = inner_marker_generator(system_params['width_local_marker']*1000, system_params['length_local_marker']*1000, layer_marker, cardinal = 'B')
        inner_marker_SE.move(destination = (np.cos(45)*system_params['length_to_local']*1000, -0.71*np.cos(45)*system_params['length_to_local']*1000))
        D_local_markers.add_ref(inner_marker_SE, alias = "inner_marker_SE")

        inner_marker_SW = inner_marker_generator(system_params['width_local_marker']*1000, system_params['length_local_marker']*1000, layer_marker, cardinal = 'A')
        inner_marker_SW.move(destination = (-np.cos(45)*system_params['length_to_local']*1000, -0.71*np.cos(45)*system_params['length_to_local']*1000))
        D_local_markers.add_ref(inner_marker_SW, alias = "inner_marker_SW")


    d_ref6 = D_local_markers.add_ref(D)
    
    """ 
    Defining the protection pad area, based on how many bonding pads there are.
    The north and south vary, but the west and east stay the same, currently.
    """
    
    prod_pads_N = [[len(bond_numbers[0]) - 1, 0]]
    prod_pads_S = [[len(bond_numbers[1]) - 1, 0]]
    prod_pads_W = [[4, 1]]
    prod_pads_E = [[4, 1]]
    
    prod_pads = [prod_pads_N, prod_pads_S, prod_pads_W, prod_pads_E]

    layer_pads = [layer_N_pads, layer_S_pads, layer_W_pads, layer_E_pads]
    
    D = pg.union(D, by_layer = True)
    
    # The dose test boolean, decides whether a dose test design to be made. The dose test
    # design is the exact same design, but with only the inner structure, and a 
    # cropped fanout. 
    
    if dose_test:
        dose_test_device = Device()
        dose_test_device.add_ref(D_without_local_markers)
        qp(dose_test_device)
        pu.write_lyp(f'2x{amount_of_dots}_dosetest.lyp', layerset = lys)

        dose_test_device.write_gds(filename = f'2x{amount_of_dots}_dosetest.gds', # Output GDS file name
                unit = 1e-9,                  # Base unit (1e-6 = microns)
                precision = 1e-9,             # Precision / resolution (1e-9 = nanometers)
                auto_rename = True,           # Automatically rename cells to avoid collisions
                max_cellname_length = 28,     # Max length of cell names
                cellname = 'toplevel'         # Name of output top-level cell
            )
    
    return D_local_markers, D_without_local_markers, list_of_numbers, ext_ports, bond_numbers, layer_pads, prod_pads


def build_and_route_bond_pads(device, list_of_numbers, bond_numbers, ext_ports, prod_pads, layer_pads, width_bonding_box = None, nested = True):
    
    """ 
    The function builds the outer structure of the chip. So all the fan-outs
    above the medium step, all the bonding pads, the protection pads, the outer makers
    the text on the chip, and the larger, visible patterns.
    
    :device: is the device
    
    :list_of_numbers: is a list of numbers deciding the placement of the bonding pads.
    
    :bond_numbers: is a list of numbers deciding the placement of the bonding pads.
    
    :ext_ports: is a list of the placement of the innermost fan-outs and their placement.
    
    :prod_pads: The protection pads meant for the bonding pads, which are not ohmic contacts.
    
    :layer_pads: A 4xX array, where X is the length of each cardinal direction of bonding pads,
    which decides the layer of the bonding pads.
    
    :width_bonding_box: The width of the bonding box. If None, the width is set to the initial
    value.
    
    :nested: boolean which decides on an angular step-like nested fan-out or straight fan-out.

    """
    
    local_marker_NE_center = pg.extract(device["inner_marker_NE0"]).polygons[0].center
    local_marker_SE_center = pg.extract(device["inner_marker_SE0"]).polygons[0].center
    local_marker_NW_center = pg.extract(device["inner_marker_NW0"]).polygons[0].center
    local_marker_SW_center = pg.extract(device["inner_marker_SW0"]).polygons[0].center
    
    
    if width_bonding_box == None:
        xmax_bonding_box = device.center[0] + system_params['width_bonding_box_init']/2*1000
        xmin_bonding_box = device.center[0] - system_params['width_bonding_box_init']/2*1000

    else:
        xmax_bonding_box = device.center[0] + width_bonding_box/2
        xmin_bonding_box = device.center[0] - width_bonding_box/2

    ymax_bonding_box = device.center[1] + system_params['height_bonding_box']/2*1000
    ymin_bonding_box = device.center[1] - system_params['height_bonding_box']/2*1000
    
    x_bond_edge = 100000
    y_bond_edge = 400000

    """ 
    Creating the arrays that will evenly distribute the bonding pads in the 
    allowed range. 
    """

    north_x0 = np.linspace(xmin_bonding_box, xmax_bonding_box - x_bond_edge, num = len(bond_numbers[0]))
    south_x0 = np.linspace(xmin_bonding_box, xmax_bonding_box - x_bond_edge, num = len(bond_numbers[1]))
    west_y0 = np.linspace(ymin_bonding_box + y_bond_edge, ymax_bonding_box - y_bond_edge, num = len(bond_numbers[2]))
    east_y0 = np.linspace(ymin_bonding_box + y_bond_edge, ymax_bonding_box - y_bond_edge, num = len(bond_numbers[3]))

    north_south_multiplier = system_params['north_south_inner_nest_width_multiplier']

    nest1_divider = 1.1
    
    nest2_divider = 5
    
    nesting_dist_1 = 1.2
    
    nesting_dist_2 = 4
    
    nest1_width = 0.7*system_params['length_bonding_pad']*1000
    
    nest2_width = 0.07*system_params['length_bonding_pad']*1000
    
    centering_rectifier = 0.5*system_params['length_bonding_pad']*1000

    nest1_north_loc = np.linspace((xmin_bonding_box + centering_rectifier)/nest1_divider, (xmax_bonding_box - centering_rectifier)/nest1_divider, num = len(bond_numbers[0]))
    nest1_south_loc = np.linspace((xmin_bonding_box + centering_rectifier)/nest1_divider, (xmax_bonding_box - centering_rectifier)/nest1_divider, num = len(bond_numbers[1]))
    nest1_west_loc = np.linspace((west_y0[0] + centering_rectifier)/nest1_divider, (west_y0[-1] + centering_rectifier)/nest1_divider, num = len(bond_numbers[2]))
    nest1_east_loc = np.linspace((east_y0[0] + centering_rectifier)/nest1_divider, (east_y0[-1] + centering_rectifier)/nest1_divider, num = len(bond_numbers[3]))

    nest2_north_loc = np.linspace((xmin_bonding_box + centering_rectifier)*north_south_multiplier, (xmax_bonding_box - centering_rectifier)*north_south_multiplier, num = len(bond_numbers[0]))
    nest2_south_loc = np.linspace((xmin_bonding_box + centering_rectifier)*north_south_multiplier, (xmax_bonding_box - centering_rectifier)*north_south_multiplier, num = len(bond_numbers[1]))
    nest2_west_loc = np.linspace(nest1_west_loc[0]/nest2_divider, nest1_west_loc[-1]/nest2_divider, num = len(bond_numbers[2]))
    nest2_east_loc = np.linspace(nest1_east_loc[0]/nest2_divider, nest1_east_loc[-1]/nest2_divider, num = len(bond_numbers[3]))
    
    
    bond_pads_N, bond_pads_S, bond_pads_W, bond_pads_E, bond_pads_port_N, bond_pads_port_S, bond_pads_port_W, bond_pads_port_E = [], [], [], [], [], [], [], []

    buffer_pad_N, buffer_pad_S, buffer_pad_W, buffer_pad_E = [], [], [], []

    nest1_N, nest1_S, nest1_W, nest1_E, nest1_port_N, nest1_port_S, nest1_port_W, nest1_port_E = [], [], [], [], [], [], [], []

    nest2_N, nest2_S, nest2_W, nest2_E, nest2_port_N, nest2_port_S, nest2_port_W, nest2_port_E = [], [], [], [], [], [], [], []

    mid_port_N, mid_port_S, mid_port_W, mid_port_E = [], [], [], []

    """ 
    This for loop is long and maybe hard to read, however, its function is
    fairly simple:
    Create the bonding pads, place them according the arrays defined above,
    create the ports in the wanted direction, and route the ports to the middle ports
    that are defined via the x_midpoint_ext_port_h and y_midpoint_ext_port_v
    functions. 
    """

    counter_1 = 0

    buffer_params = parameters["buffer-params"]                     # Perhaps this is an outdates parameter, as the new radius is simply slightly smaller than the bonding pad width
    
    buffer_radius = system_params['width_bonding_pad']*1000/2*0.9

    buffer_layer = metal_buffer_layer[1]

    east_west_inner_nest_width_multiplier = system_params['east_west_inner_nest_width_multiplier']

    distance_ext_port_north = system_params['distance_external_port_north']*1000
    distance_ext_port_south = system_params['distance_external_port_south']*1000
    distance_ext_port_west = system_params['distance_external_port_west']*1000
    distance_ext_port_east = system_params['distance_external_port_east']*1000

    distance_nest_2_north = system_params['distance_nest_2_north']*1000
    distance_nest_2_south = system_params['distance_nest_2_south']*1000
    distance_nest_2_west = system_params['distance_nest_2_west']*1000
    distance_nest_2_east = system_params['distance_nest_2_east']*1000

    xy_spots_multiplier = 200

    for i in list_of_numbers:
        for j in range(len(i)):
            if len(i) % 2 == 0:
                x_y_spots_init = np.arange(-len(i)/2, len(i)/2 + 1, 1 , dtype = int)
                x_y_spots = np.delete(x_y_spots_init, np.where(x_y_spots_init == 0))
            else:
                x_y_spots= np.arange(-math.floor(len(i)/2), math.floor(len(i)/2) + 1, 1 , dtype = int)
            if counter_1 == 0:
                bond_pads_N.append("bonding_pad_N" + str(i[j]))
                bond_pads_port_N.append("bond_port_pad_N" + str(i[j]))
                bond_pads_N[j] = device << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[0][j]).rotate(90)
                bond_pads_N[j].ymin = ymax_bonding_box
                bond_pads_N[j].xmin = north_x0[j]
                buffer_pad_N.append("buffer_pad_N" + str(i[j]))
                buffer_pad_N[j] = device << pg.circle(radius = buffer_radius, layer = buffer_layer)
                buffer_pad_N[j].center = [north_x0[j] + 0.5*width_bonding_pad, bond_pads_N[j].ymax - 2*0.55*buffer_radius]
                bond_pads_port_N[j] = device.add_port(name='B_N' + str(i[j]), midpoint = [bond_pads_N[j].center[0], bond_pads_N[j].ymin], width = system_params['width_bonding_pad']*1000, orientation = 270)
                if nested:
                    nest1_N.append("nest1_N" + str(i[j]))
                    nest1_port_N.append("nest1_port_N" + str(i[j]))
                    nest2_N.append("nest2_N" + str(i[j]))
                    nest2_port_N.append("nest2_port_N" + str(i[j]))
                    nest1_port_N[j] = device.add_port(name='nest1_N' + str(i[j]), midpoint = [nest1_north_loc[j], ymax_bonding_box/nesting_dist_1], width = nest1_width, orientation = 270)
                    nest2_port_N[j] = device.add_port(name='nest2_N' + str(i[j]), midpoint = [nest2_north_loc[j], distance_nest_2_north], width = nest2_width, orientation = 270)
                    device.add_ref(pr.route_quad(bond_pads_port_N[j], nest1_port_N[j], layer = layer_pads[0][j]))
                    device.add_ref(pr.route_quad(nest1_port_N[j], nest2_port_N[j], layer = layer_pads[0][j]))
                    device.add_ref(pr.route_quad(nest2_port_N[j], ext_ports[0][j], layer = layer_pads[0][j]))
                else:
                    device.add_ref(pr.route_quad(bond_pads_port_N[j], ext_ports[0][j], layer = layer_pads[0][j]))
                mid_port_N.append('mid_port_N' + str(j))
                mid_port_N[j] = device.add_port(name='MID_N' + str(j),
                            midpoint = [x_midpoint_ext_port_h(i[j]) - xy_spots_multiplier*x_y_spots[j], distance_ext_port_north - system_params['closing_distance']*1000],
                            width = system_params['width_middle_port']*1000, orientation = 270)
                device.add_ref(pr.route_quad(ext_ports[0][j], mid_port_N[j], layer = layer_pads[0][j]))
            if counter_1 == 1:
                bond_pads_S.append("bonding_pad_S" + str(i[j]))
                bond_pads_port_S.append("bond_port_pad_S" + str(i[j]))
                bond_pads_S[j] = device << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[1][j]).rotate(90)
                bond_pads_S[j].ymax = ymin_bonding_box
                bond_pads_S[j].xmin = south_x0[j]
                buffer_pad_S.append("buffer_pad_S" + str(i[j]))
                buffer_pad_S[j] = device << pg.circle(radius = buffer_radius, layer = buffer_layer)
                buffer_pad_S[j].center = [south_x0[j] + 0.5*width_bonding_pad, bond_pads_S[j].ymin + 2*0.55*buffer_radius]
                bond_pads_port_S[j] = device.add_port(name='B_S' + str(i[j]),
                            midpoint = [bond_pads_S[j].center[0], bond_pads_S[j].ymax], width = system_params['width_bonding_pad']*1000, orientation = 90)
                if nested:
                    nest1_S.append("nest1_S" + str(i[j]))
                    nest1_port_S.append("nest1_port_S" + str(i[j]))
                    nest2_S.append("nest2_S" + str(i[j]))
                    nest2_port_S.append("nest2_port_S" + str(i[j]))
                    nest1_port_S[j] = device.add_port(name='nest1_S' + str(i[j]), midpoint = [nest1_south_loc[j], ymin_bonding_box/nesting_dist_1], width = nest1_width, orientation = 90)
                    nest2_port_S[j] = device.add_port(name='nest2_S' + str(i[j]), midpoint = [nest2_south_loc[j], distance_nest_2_south], width = nest2_width, orientation = 90)
                    device.add_ref(pr.route_quad(bond_pads_port_S[j], nest1_port_S[j], layer = layer_pads[1][j]))
                    device.add_ref(pr.route_quad(nest1_port_S[j], nest2_port_S[j], layer = layer_pads[1][j]))
                    device.add_ref(pr.route_quad(nest2_port_S[j], ext_ports[1][j], layer = layer_pads[1][j]))
                else:
                    device.add_ref(pr.route_quad(bond_pads_port_S[j], ext_ports[1][j], layer = layer_pads[1][j]))
                mid_port_S.append('mid_port_S' + str(j))
                mid_port_S[j] = device.add_port(name='MID_S' + str(j),
                            midpoint = [x_midpoint_ext_port_h(i[j]) - xy_spots_multiplier*x_y_spots[j], distance_ext_port_south + system_params['closing_distance']*1000],
                            width = system_params['width_middle_port']*1000, orientation = 90)
                device.add_ref(pr.route_quad(ext_ports[1][j], mid_port_S[j], layer = layer_pads[1][j]))
            if counter_1 == 2:
                bond_pads_W.append("bonding_pad_W" + str(i[j]))
                bond_pads_port_W.append("bond_port_pad_W" + str(i[j]))
                bond_pads_W[j] = device << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[2][j])
                bond_pads_W[j].xmax = xmin_bonding_box
                bond_pads_W[j].ymin = west_y0[j]
                if j > 0 and j < len(i) - 1:
                    buffer_pad_W.append("buffer_pad_W" + str(i[j]))
                    buffer_pad_W[j - 1] = device << pg.circle(radius = buffer_radius, layer = buffer_layer)
                    buffer_pad_W[j - 1].center = [bond_pads_W[j].xmin + 2*0.55*buffer_radius, west_y0[j] + 0.5*width_bonding_pad]
                bond_pads_port_W[j] = device.add_port(name='B_W' + str(i[j]), midpoint = [bond_pads_W[j].xmax, bond_pads_W[j].center[1]],
                            width = system_params['width_bonding_pad']*1000, orientation = 0)
                if nested:
                    nest1_W.append("nest1_W" + str(i[j]))
                    nest1_port_W.append("nest1_port_W" + str(i[j]))
                    nest2_W.append("nest2_W" + str(i[j]))
                    nest2_port_W.append("nest2_port_W" + str(i[j]))
                    nest1_port_W[j] = device.add_port(name='nest1_W' + str(i[j]), midpoint = [xmin_bonding_box/nesting_dist_1, nest1_west_loc[j]], width = nest1_width, orientation = 0)
                    nest2_port_W[j] = device.add_port(name='nest2_W' + str(i[j]), midpoint = [distance_nest_2_west, nest2_west_loc[j]], width = nest2_width, orientation = 0)
                    device.add_ref(pr.route_quad(bond_pads_port_W[j], nest1_port_W[j], layer = layer_pads[2][j]))
                    device.add_ref(pr.route_quad(nest1_port_W[j], nest2_port_W[j], layer = layer_pads[2][j]))
                    device.add_ref(pr.route_quad(nest2_port_W[j], ext_ports[2][j], layer = layer_pads[2][j]))
                else:
                    device.add_ref(pr.route_quad(bond_pads_port_W[j], ext_ports[2][j], layer = layer_pads[2][j]))
                mid_port_W.append('mid_port_W' + str(j))
                mid_port_W[j] = device.add_port(name='MID_W' + str(j),
                            midpoint = [distance_ext_port_west + system_params['closing_distance']*1000, y_midpoint_ext_port_v(i[j]) - xy_spots_multiplier*x_y_spots[j]],
                            width = system_params['width_middle_port']*1000, orientation = 0)
                device.add_ref(pr.route_quad(ext_ports[2][j], mid_port_W[j], layer = layer_pads[2][j]))
            if counter_1 == 3:
                bond_pads_E.append("bonding_pad_E" + str(i[j]))
                bond_pads_port_E.append("bond_port_pad_E" + str(i[j]))
                bond_pads_E[j] = device << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[3][j])
                bond_pads_E[j].xmin = xmax_bonding_box
                bond_pads_E[j].ymin = east_y0[j]
                if j > 0 and j < len(i) - 1:
                    buffer_pad_E.append("buffer_pad_E" + str(i[j]))
                    buffer_pad_E[j - 1] = device << pg.circle(radius = buffer_radius, layer = buffer_layer)
                    buffer_pad_E[j - 1].center = [bond_pads_E[j].xmax - 2*0.55*buffer_radius, east_y0[j] + 0.5*width_bonding_pad]
                bond_pads_port_E[j] = device.add_port(name='B_E' + str(i[j]), midpoint = [bond_pads_E[j].xmin, bond_pads_E[j].center[1]],
                            width = system_params['width_bonding_pad']*1000, orientation = 180)
                if nested:
                    nest1_E.append("nest1_E" + str(i[j]))
                    nest1_port_E.append("nest1_port_E" + str(i[j]))
                    nest2_E.append("nest2_E" + str(i[j]))
                    nest2_port_E.append("nest2_port_E" + str(i[j]))
                    nest1_port_E[j] = device.add_port(name='nest1_E' + str(i[j]), midpoint = [xmax_bonding_box/nesting_dist_1, nest1_east_loc[j]], width = nest1_width, orientation = 180)
                    nest2_port_E[j] = device.add_port(name='nest2_E' + str(i[j]), midpoint = [distance_nest_2_east, nest2_east_loc[j]], width = nest2_width, orientation = 180)
                    device.add_ref(pr.route_quad(bond_pads_port_E[j], nest1_port_E[j], layer = layer_pads[3][j]))
                    device.add_ref(pr.route_quad(nest1_port_E[j], nest2_port_E[j], layer = layer_pads[3][j]))
                    device.add_ref(pr.route_quad(nest2_port_E[j], ext_ports[3][j], layer = layer_pads[3][j]))
                else:
                    device.add_ref(pr.route_quad(bond_pads_port_E[j], ext_ports[3][j], layer = layer_pads[3][j]))
                mid_port_E.append('mid_port_E' + str(j))
                mid_port_E[j] = device.add_port(name='MID_E' + str(j),
                            midpoint = [distance_ext_port_east - system_params['closing_distance']*1000, y_midpoint_ext_port_v(i[j]) - xy_spots_multiplier*x_y_spots[j]],
                            width = system_params['width_middle_port']*1000, orientation = 180)
                device.add_ref(pr.route_quad(ext_ports[3][j], mid_port_E[j], layer = layer_pads[3][j]))
        counter_1 += 1

    counter_1 = 0

    prod_pads_bool = True

    """ 
    Given the boolean prod_pads_bool, create protection pads for the bonding
    pads in each cardinal direction. 
    """

    if prod_pads_bool:
        for i in prod_pads:
            for j in range(len(i)):
                if counter_1 == 0:
                    protection_pads_N = protection_pads(device, (bond_pads_N[i[j][0]].xmax + protection_extra - (bond_pads_N[i[j][1]].xmin - protection_extra)), system_params['length_bonding_pad']*1000,
                        destination = (bond_pads_N[i[j][1]].xmin - protection_extra,
                                       bond_pads_N[i[j][1]].ymin + protection_extra))
                if counter_1 == 1:
                    protection_pads_S = protection_pads(device, (bond_pads_S[i[j][0]].xmax + protection_extra - (bond_pads_S[i[j][1]].xmin - protection_extra)), system_params['length_bonding_pad']*1000,
                        destination = (bond_pads_S[i[j][1]].xmin - protection_extra,
                                       bond_pads_S[i[j][1]].ymin - protection_extra))
                if counter_1 == 2:
                    protection_pads_W = protection_pads(device, system_params['length_bonding_pad']*1000, (bond_pads_W[i[j][0]].ymax + protection_extra - (bond_pads_W[i[j][1]].ymin - protection_extra)),
                        destination = (bond_pads_W[i[j][1]].xmin - protection_extra,
                                       bond_pads_W[i[j][1]].ymin - protection_extra))
                if counter_1 == 3:
                    protection_pad_E = protection_pads(device, system_params['length_bonding_pad']*1000, (bond_pads_E[i[j][0]].ymax + protection_extra - (bond_pads_E[i[j][1]].ymin - protection_extra)),
                        destination = (bond_pads_E[i[j][1]].xmin + protection_extra,
                                       bond_pads_E[i[j][1]].ymin - protection_extra))
            counter_1 += 1

    device_joined_by_layer = pg.union(device, by_layer = True)

    return device_joined_by_layer


def join_structures_on_chip(devices, sem_devices, layer_dummy_bonds, full = True):
    
    """ 
    Make the full device by joining the structures on the chip. Then creates 
    all the various other things, such as the dummy pads, the visual markers, 
    the name on the chip.
    
    :outer_params: The set of parameters, which defines how the chip will look. If None,
    then the default parameters, defined above, will be used.
    
    :devices: is a list of the devices you want on the chip. It is set to be 4 long.
    
    :sem_devices: is a list of the SEM structures you want on the chip. It is set to be 4 long.
    
    :bottom_left: is a list of the relative coordinates for the bottom left 
    grouping of markers. It is defined such that the corner marker is at position 0,0
    and the rest are at integer position in relation to that.
    
    :bottom_right: is a list of the relative coordinates for the bottom right 
    grouping of markers. It is defined such that the corner marker is at position 0,0
    and the rest are at integer position in relation to that.
    
    :top_left: read above
    
    :top_right_ read above.
    """
    
    outer_params = parameters['outer-params']
    distance_text_from_border = outer_params['distance_text_from_y_border']*1000
    distance_text_from_middle = outer_params['distance_text_from_middle']*1000
    visual_marker_length = outer_params['visual_marker_length']*1000
    visual_marker_width = outer_params['visual_marker_width']*1000
    distance_structure_from_center = outer_params['distance_structure_from_center']
    width_die = outer_params['width_die']*1000000
    length_die = outer_params['length_die']*1000000
    name_of_chip = outer_params['name_of_chip']
    layer_slice_visuals = layer_marker
    x_top_left = marker_params['x_top_left']
    y_top_left = marker_params['y_top_left']
    x_top_right = marker_params['x_top_right']
    y_top_right = marker_params['y_top_right']
    x_bottom_left = marker_params['x_bottom_left']
    y_bottom_left = marker_params['y_bottom_left']
    x_bottom_right = marker_params['x_bottom_right']
    y_bottom_right = marker_params['y_bottom_right']
    distance_from_border_marker = marker_params['distance_from_border']*1000
    
    top_left = [x_top_left, y_top_left]
    top_right = [x_top_right, y_top_right]
    bottom_left = [x_bottom_left, y_bottom_left]
    bottom_right = [x_bottom_right, y_bottom_right]
    
    if not full:
        length_die = 3.2e6
    
    if outer_params == None:
        outer_params = outer_params
    full_device = Device()
    
    full_refs = np.zeros(4, dtype = tuple)
    
    die = pg.rectangle(size = (length_die, width_die), layer = layer_die)
    
    """ 
    Place all the structures on the device, subsequently move their positions
    and create and place the name of the chip. 
    """
    
    names = ["A0", "B0", "C0", "D0", "A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2", "A3", "B3", "C3", "D3"]
    
    for i in range(len(devices)):
        full_refs[i] = full_device.add_ref(devices[i])


    full_refs[0].xmin = die.center[0] - distance_structure_from_center*width_die - 0.5*(full_refs[0].xmax - full_refs[0].xmin)
    full_refs[0].ymax = die.center[1] - distance_structure_from_center*width_die + 0.5*(full_refs[0].ymax - full_refs[0].ymin)        
            
    full_refs[1].xmax = die.center[0] - distance_structure_from_center*width_die + 0.5*(full_refs[1].xmax - full_refs[1].xmin)
    full_refs[1].ymin = die.center[1] + distance_structure_from_center*width_die - 0.5*(full_refs[1].ymax - full_refs[1].ymin)


    if full:
        full_refs[2].xmin = die.center[0] + distance_structure_from_center*width_die - 0.5*(full_refs[2].xmax - full_refs[2].xmin)
        full_refs[2].ymin = die.center[1] + distance_structure_from_center*width_die - 0.5*(full_refs[2].ymax - full_refs[2].ymin)
        
        full_refs[3].xmin = die.center[0] + distance_structure_from_center*width_die - 0.5*(full_refs[3].xmax - full_refs[3].xmin)
        full_refs[3].ymax = die.center[1] - distance_structure_from_center*width_die + 0.5*(full_refs[3].ymax - full_refs[3].ymin)

    if not full:
        full_refs[0].xmin += 0.5*(full_refs[0].xmax - full_refs[0].xmin)*1.1
        full_refs[1].xmax += 0.5*(full_refs[1].xmax - full_refs[1].xmin)*1.1
    

    transform = system_params['transform']

    local_marker_positions = {}

    local_marker_pos_1 = inner_params['local_marker_positions_1']
    local_marker_pos_2 = inner_params['local_marker_positions_3']
    local_marker_pos_3 = inner_params['local_marker_positions_2']
    local_marker_pos_4 = inner_params['local_marker_positions_4']
    
    local_marker_pos = [local_marker_pos_1, local_marker_pos_2, local_marker_pos_3, local_marker_pos_4]

    correct_order_refs = [full_refs[1], full_refs[2], full_refs[0], full_refs[3]]
    
    k = 0
    for i in range(len(full_refs)):
        for j in range(4):
            local_positions = [[round((correct_order_refs[i].center[0] - math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12), round((correct_order_refs[i].center[1] - math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12)],
                         [round((correct_order_refs[i].center[0] + math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12), round((correct_order_refs[i].center[1] - math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12)],
                         [round((correct_order_refs[i].center[0] - math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12), round((correct_order_refs[i].center[1] + math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12)],
                         [round((correct_order_refs[i].center[0] + math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12), round((correct_order_refs[i].center[1] + math.sqrt(1/2)*local_marker_pos[i][j]*1000)/1000000 + transform, 12)]]
            local_marker_positions.update({str(i + 1) + names[k]: local_positions[0]})
            k += 1
            local_marker_positions.update({str(i + 1) + names[k]: local_positions[1]})
            k += 1
            local_marker_positions.update({str(i + 1) + names[k]: local_positions[2]})
            k += 1
            local_marker_positions.update({str(i + 1) + names[k]: local_positions[3]})
            k += 1
            write_SCON("LM_" + str(i + 1) + str(j), local_positions, "local")
        k = 0

    if not os.path.isdir('Marker Positions'):
        os.mkdir('Marker Positions')
    
    file_name_for_local_markers = working_directory + "/Marker Positions/" + "Local_marker_placements_" + name_of_chip + ".json"
    with open(file_name_for_local_markers, 'w') as f:
        json.dump(local_marker_positions, f)
    
    text_die = full_device << pg.text(text=name_of_chip, size=1.8e5, justify='left', layer=layer_text, font='DEPLOF')
    text_die.ymax = die.ymax - distance_text_from_border
    text_die.xmin = die.center[0] + distance_text_from_middle

    if not full:
        text_die.xmin = die.center[0] - distance_text_from_border 
        

    """ 
    Generate and place the markers in each corner of the chip.
    """
    
    length_outer_marker = system_params['length_outer_marker']*1000
    length_inner_marker = system_params['length_inner_marker']*1000
    
    top_left_markers_device, top_left_markers, top_left_group, top_left_names = marker_device_generator(full_device, top_left[0], top_left[1], 'C')
    top_right_markers_device, top_right_markers, top_right_group, top_right_names = marker_device_generator(full_device, top_right[0], top_right[1], 'D')
    bottom_left_markers_device, bottom_left_markers, bottom_left_group, bottom_left_names = marker_device_generator(full_device, bottom_left[0], bottom_left[1], 'A')
    bottom_right_markers_device, bottom_right_markers, bottom_right_group, bottom_right_names = marker_device_generator(full_device, bottom_right[0], bottom_right[1], 'B')

    top_left_group.xmin = die.xmin + distance_from_border_marker + 0.5*length_outer_marker - 0.5*length_inner_marker
    top_left_group.ymax = die.ymax - distance_from_border_marker - 0.5*length_outer_marker + 0.5*length_inner_marker

    top_right_group.xmax = die.xmax - distance_from_border_marker - 0.5*length_outer_marker + 0.5*length_inner_marker
    top_right_group.ymax = die.ymax - distance_from_border_marker - 0.5*length_outer_marker + 0.5*length_inner_marker

    bottom_right_group.xmax = die.xmax - distance_from_border_marker - 0.5*length_outer_marker + 0.5*length_inner_marker
    bottom_right_group.ymin = die.ymin + distance_from_border_marker + 0.5*length_outer_marker - 0.5*length_inner_marker

    bottom_left_group.xmin = die.xmin + distance_from_border_marker + 0.5*length_outer_marker - 0.5*length_inner_marker
    bottom_left_group.ymin = die.ymin + distance_from_border_marker + 0.5*length_outer_marker - 0.5*length_inner_marker
    
    top_left_markers_device.xmin = die.xmin + distance_from_border_marker
    top_left_markers_device.ymax = die.ymax - distance_from_border_marker

    top_right_markers_device.xmax = die.xmax - distance_from_border_marker 
    top_right_markers_device.ymax = die.ymax - distance_from_border_marker

    bottom_right_markers_device.xmax = die.xmax - distance_from_border_marker 
    bottom_right_markers_device.ymin = die.ymin + distance_from_border_marker

    bottom_left_markers_device.xmin = die.xmin + distance_from_border_marker
    bottom_left_markers_device.ymin = die.ymin + distance_from_border_marker
    
    for i in range(len(top_left_markers)):
        full_device.add_ref(top_left_markers[i])
        full_device.add_ref(top_right_markers[i])
        full_device.add_ref(bottom_left_markers[i])
        full_device.add_ref(bottom_right_markers[i])

    full_device.add_ref(top_left_markers_device)
    full_device.add_ref(top_right_markers_device)
    full_device.add_ref(bottom_left_markers_device)
    full_device.add_ref(bottom_right_markers_device)

    global_marker_positions = {}
    
    for i in range(len(top_left_names)):
        global_positions = [[bottom_left_markers[i].center[0]/1000000 + transform, bottom_left_markers[i].center[1]/1000000 + transform],
                            [bottom_right_markers[i].center[0]/1000000 + transform, bottom_right_markers[i].center[1]/1000000 + transform],
                            [top_left_markers[i].center[0]/1000000 + transform, top_left_markers[i].center[1]/1000000 + transform],
                            [top_right_markers[i].center[0]/1000000 + transform, top_right_markers[i].center[1]/1000000 + transform]]
        
        global_marker_positions.update({bottom_left_names[i]: global_positions[0]})
        global_marker_positions.update({bottom_right_names[i]: global_positions[1]})
        global_marker_positions.update({top_left_names[i]: global_positions[2]})
        global_marker_positions.update({top_right_names[i]: global_positions[3]})
        
        write_SCON("GM_" + bottom_left_names[i] + bottom_right_names[i] + top_left_names[i] + top_right_names[i], global_positions, "global")
    
    file_name_for_global_markers = working_directory + "/Marker Positions/" + "Global_marker_placements_" + name_of_chip + ".json"
    with open(file_name_for_global_markers, 'w') as f:
        json.dump(global_marker_positions, f)

    """ 
    Set parameters for the dummy bond pads and utilise the dummy_pad function.
    
    See above documentation for further information on the 
    :create_and_place_dumme_bond_pads: function.
    """
    
    spacing_dummy_bonding_pads = 130e3 #nm

    amount_dummy_bonding_pads = 3

    dummy_pads_distance_from_center = 250e3

    dummy_pads_distance_from_side = 200e3
    
    if full:
        dummy_pads_NN = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads, 
                length_bonding_pad, width_bonding_pad, spacing_dummy_bonding_pads, 'x',
            destination = (die.center[0] - 0.1*dummy_pads_distance_from_center - 5*spacing_dummy_bonding_pads - 6*width_bonding_pad, die.ymax - 2*dummy_pads_distance_from_center), layer_dummy_bonds = layer_dummy_bonds)
        dummy_pads_SSW = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
                length_bonding_pad, width_bonding_pad, spacing_dummy_bonding_pads, 'x',
            destination = (die.center[0] - 0.1*dummy_pads_distance_from_center - 2*spacing_dummy_bonding_pads - 6*width_bonding_pad, die.ymin + 1.5*dummy_pads_distance_from_center), layer_dummy_bonds = layer_dummy_bonds)
        dummy_pads_SSE = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
                length_bonding_pad, width_bonding_pad, spacing_dummy_bonding_pads, 'x',
            destination = (die.center[0] + 1.5*dummy_pads_distance_from_center, die.ymin + 1.5*dummy_pads_distance_from_center), layer_dummy_bonds = layer_dummy_bonds)
    dummy_pads_SE = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
            width_bonding_pad, length_bonding_pad, spacing_dummy_bonding_pads, 'y',
            destination = (die.xmax - 2*dummy_pads_distance_from_side - length_bonding_pad, die.center[1] - 0.1*dummy_pads_distance_from_center - 2*spacing_dummy_bonding_pads - 6*width_bonding_pad), layer_dummy_bonds = layer_dummy_bonds)
    dummy_pads_NE = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
            width_bonding_pad, length_bonding_pad, spacing_dummy_bonding_pads, 'y',
            destination = (die.xmax - 2*dummy_pads_distance_from_side  - length_bonding_pad, die.center[1] + dummy_pads_distance_from_center), layer_dummy_bonds = layer_dummy_bonds)
    dummy_pads_NW = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
            width_bonding_pad, length_bonding_pad, spacing_dummy_bonding_pads, 'y',
            destination = ((die.xmin + 2*dummy_pads_distance_from_side, die.center[1] + dummy_pads_distance_from_center)), layer_dummy_bonds = layer_dummy_bonds)
    dummy_pads_SW = create_and_place_dummy_bond_pads(full_device, amount_dummy_bonding_pads,
            width_bonding_pad, length_bonding_pad, spacing_dummy_bonding_pads, 'y',
            destination = (die.xmin + 2*dummy_pads_distance_from_side, die.center[1] - 0.1*dummy_pads_distance_from_center - 2*spacing_dummy_bonding_pads - 6*width_bonding_pad), layer_dummy_bonds = layer_dummy_bonds)
    """ 
    Define the bond pads that need protection, iterate through all cardinal directions,
    and create a protectional layer for the barriers, dots, and screening gates.
    The ohmics are not given protectional layers, as it does not matter.
    See above documentation for more information about the :protection_pads: function.
    """
    if full:
        for_protection_N = [dummy_pads_NN]

        for_protection_SS = [dummy_pads_SSW, dummy_pads_SSE]

    for_protection_W = [dummy_pads_SW, dummy_pads_NW]

    for_protection_E = [dummy_pads_SE, dummy_pads_NE]
        
    if full: 
        for pads in for_protection_N:
            protection_pad = protection_pads(full_device, (pads[2].xmax + protection_extra - (pads[1].xmin - protection_extra)), length_bonding_pad,
                    destination = (pads[1].xmin - protection_extra, pads[1].ymin + protection_extra))
            # protection_pad = protection_pads(full_device, (pads[5].xmax + protection_extra - (pads[4].xmin - protection_extra)), length_bonding_pad,
            #         destination = (pads[4].xmin - protection_extra, pads[4].ymin + protection_extra))
            
        for pads in for_protection_SS:
            protection_pad = protection_pads(full_device, (pads[2].xmax + protection_extra - (pads[1].xmin - protection_extra)), length_bonding_pad,
                    destination = (pads[1].xmin - protection_extra, pads[1].ymin - protection_extra))
            # protection_pad = protection_pads(full_device, (pads[5].xmax + protection_extra - (pads[4].xmin - protection_extra)), length_bonding_pad,
            #         destination = (pads[4].xmin - protection_extra, pads[4].ymin - protection_extra))
        
    for pads in for_protection_E:
        protection_pad = protection_pads(full_device, (pads[2].xmax - protection_extra - (pads[1].xmin - protection_extra)), 2*width_bonding_pad + spacing_dummy_bonding_pads + 2*protection_extra,
                    destination = (pads[1].xmin + protection_extra, pads[1].ymin - protection_extra))
        # protection_pad = protection_pads(full_device, (pads[5].xmax - protection_extra - (pads[4].xmin - protection_extra)), 2*width_bonding_pad + spacing_dummy_bonding_pads + 2*protection_extra,
        #             destination = (pads[4].xmin + protection_extra, pads[4].ymin - protection_extra))
            

    for pads in for_protection_W:
        protection_pad = protection_pads(full_device, (pads[2].xmax - protection_extra - (pads[1].xmin - protection_extra)), 2*width_bonding_pad + spacing_dummy_bonding_pads + 2*protection_extra,
                    destination = (pads[1].xmin - protection_extra, pads[1].ymin - protection_extra))
        # protection_pad = protection_pads(full_device, (pads[5].xmax - protection_extra - (pads[4].xmin - protection_extra)), 2*width_bonding_pad + spacing_dummy_bonding_pads + 2*protection_extra,
        #             destination = (pads[4].xmin - protection_extra, pads[4].ymin - protection_extra))


    transform = system_params['transform']

    """ 
    Adding the SEM structures. Their x- and y-placement is set manually.
    """
    if full:
        sw_marker1, se_marker1, nw_marker1, ne_marker1 = place_sem_structure(full_device, sem_devices[1], 'SW', xmin = die.center[0] - 0.5*sem_devices[0].xsize, ymin =  die.center[1] + 1.8*sem_devices[0].ysize)
        sw_marker2, se_marker2, nw_marker2, ne_marker2 = place_sem_structure(full_device, sem_devices[2], 'SW', xmin = die.center[0] - 0.5*sem_devices[1].xsize, ymin =  die.center[1] + 0.5*sem_devices[1].ysize)
        sw_marker3, se_marker3, nw_marker3, ne_marker3 = place_sem_structure(full_device, sem_devices[0], 'SW', xmin = die.center[0] - 0.5*sem_devices[2].xsize, ymin =  die.center[1] - 0.8*sem_devices[2].ysize)
        sw_marker4, se_marker4, nw_marker4, ne_marker4 = place_sem_structure(full_device, sem_devices[3], 'SW', xmin = die.center[0] - 0.5*sem_devices[3].xsize, ymin =  die.center[1] - 2.1*sem_devices[3].ysize)
    
        sem_positions_1 = [[sw_marker1.center[0]/1000000 + transform, sw_marker1.center[1]/1000000 + transform],
                            [se_marker1.center[0]/1000000 + transform, se_marker1.center[1]/1000000 + transform],
                            [nw_marker1.center[0]/1000000 + transform, nw_marker1.center[1]/1000000 + transform],
                            [ne_marker1.center[0]/1000000 + transform, ne_marker1.center[1]/1000000 + transform]]
        sem_positions_2 =   [[sw_marker2.center[0]/1000000 + transform, sw_marker2.center[1]/1000000 + transform],
                            [se_marker2.center[0]/1000000 + transform, se_marker2.center[1]/1000000 + transform],
                            [nw_marker2.center[0]/1000000 + transform, nw_marker2.center[1]/1000000 + transform],
                            [ne_marker2.center[0]/1000000 + transform, ne_marker2.center[1]/1000000 + transform]]
        sem_positions_3 =   [[sw_marker3.center[0]/1000000 + transform, sw_marker3.center[1]/1000000 + transform],
                            [se_marker3.center[0]/1000000 + transform, se_marker3.center[1]/1000000 + transform],
                            [nw_marker3.center[0]/1000000 + transform, nw_marker3.center[1]/1000000 + transform],
                            [ne_marker3.center[0]/1000000 + transform, ne_marker3.center[1]/1000000 + transform]]
        sem_positions_4 =   [[sw_marker4.center[0]/1000000 + transform, sw_marker4.center[1]/1000000 + transform],
                            [se_marker4.center[0]/1000000 + transform, se_marker4.center[1]/1000000 + transform],
                            [nw_marker4.center[0]/1000000 + transform, nw_marker4.center[1]/1000000 + transform],
                            [ne_marker4.center[0]/1000000 + transform, ne_marker4.center[1]/1000000 + transform]]
        
        sem_positions = [sem_positions_1, sem_positions_2, sem_positions_3, sem_positions_4]
    
        sem_marker_placements = {
            'SW1': sem_positions_1[0],
            'SE1': sem_positions_1[1],
            'NW1': sem_positions_1[2],
            'NE1': sem_positions_1[3],
            'SW2': sem_positions_2[0],
            'SE2': sem_positions_2[1],
            'NW2': sem_positions_2[2],
            'NE2': sem_positions_2[3],
            'SW3': sem_positions_3[0],
            'SE3': sem_positions_3[1],
            'NW3': sem_positions_3[2],
            'NE3': sem_positions_3[3],
            'SW4': sem_positions_4[0],
            'SE4': sem_positions_4[1],
            'NW4': sem_positions_4[2],
            'NE4': sem_positions_4[3],
        }
        
        for coords in sem_positions:
            write_SCON("SEM_" + str(sem_positions.index(coords)), coords, "local")
            
            
        
        file_name_for_sem_markers = working_directory + "/Marker Positions/" + "SEM_marker_placements_" + name_of_chip + ".json"
        with open(file_name_for_sem_markers, 'w') as f:
            json.dump(sem_marker_placements, f)
    """ 
    Adding visual markers.
    """

    visual_marker_west = full_device << pg.rectangle(size = (visual_marker_length, visual_marker_width), layer = layer_slice_visuals)
    visual_marker_west.move(destination = (die.center[0] - 0.43*width_die, die.center[1] - 0.5*visual_marker_width))

    visual_marker_east = full_device << pg.rectangle(size = (visual_marker_length, visual_marker_width), layer = layer_slice_visuals)
    visual_marker_east.move(destination = (die.center[0] + 0.43*width_die - visual_marker_length, die.center[1] - 0.5*visual_marker_width))

    if not full:
        visual_marker_west.move(destination = (die.center[0], die.center[1] - 0.5*width_die))
        visual_marker_east.move(destination = (die.center[0] - 0.53*width_die, die.center[1] - 0.5*width_die))
    
    if full:

        visual_marker_south = full_device << pg.rectangle(size = (visual_marker_width, visual_marker_length), layer = layer_slice_visuals)
        visual_marker_south.move(destination = (die.center[0] - 0.5*visual_marker_width, die.center[1] - 0.43*length_die))

        visual_marker_north = full_device << pg.rectangle(size = (visual_marker_width, 0.8*visual_marker_length), layer = layer_slice_visuals)
        visual_marker_north.move(destination = (die.center[0] - 0.50*visual_marker_width, die.center[1] + 0.43*length_die - 1.5*visual_marker_length))

    full_device = pg.union(full_device, by_layer = True)
    
    """
    Saving the design, along with it layer properties, and the parameters that 
    created it. We save the parameters using pickling.
    """

    if not os.path.isdir('Designs'):
        os.mkdir('Designs')
    
    pu.write_lyp(f'Designs/{name_of_chip}_layerinfo.lyp', layerset = lys)


    full_device.write_gds(filename = f'Designs/{name_of_chip}.gds', # Output GDS file name
                unit = 1e-9,                  # Base unit (1e-6 = microns)
                precision = 1e-9,             # Precision / resolution (1e-9 = nanometers)
                auto_rename = True,           # Automatically rename cells to avoid collisions
                max_cellname_length = 28,     # Max length of cell names
                cellname = 'toplevel'         # Name of output top-level cell
            )


def join_structures_test(devices, full = True):
    
    """ 
    Make the full device by joining the structures on the chip. Then creates 
    all the various other things, such as the dummy pads, the visual markers, 
    the name on the chip.
    
    :outer_params: The set of parameters, which defines how the chip will look. If None,
    then the default parameters, defined above, will be used.
    
    :devices: is a list of the devices you want on the chip. It is set to be 4 long.
    
    :sem_devices: is a list of the SEM structures you want on the chip. It is set to be 4 long.
    
    :bottom_left: is a list of the relative coordinates for the bottom left 
    grouping of markers. It is defined such that the corner marker is at position 0,0
    and the rest are at integer position in relation to that.
    
    :bottom_right: is a list of the relative coordinates for the bottom right 
    grouping of markers. It is defined such that the corner marker is at position 0,0
    and the rest are at integer position in relation to that.
    
    :top_left: read above
    
    :top_right_ read above.
    """
    
    outer_params = parameters['outer-params']
    distance_text_from_border = outer_params['distance_text_from_y_border']*1000
    distance_text_from_middle = outer_params['distance_text_from_middle']*1000
    distance_structure_from_center = outer_params['distance_structure_from_center']
    width_die = outer_params['width_die']*1000000
    length_die = outer_params['length_die']*1000000
    name_of_chip = outer_params['name_of_chip']

    
    if not full:
        length_die = 3.2e6
    
    if outer_params == None:
        outer_params = outer_params
    full_device = Device()
    
    full_refs = np.zeros(4, dtype = tuple)
    
    die = pg.rectangle(size = (length_die, width_die), layer = layer_die)
    
    """ 
    Place all the structures on the device, subsequently moce their positions
    and create and place the name of the chip. 
    """
    
    full_refs = full_device.add_ref(devices[0])

    full_refs.xmin = die.center[0] - distance_structure_from_center*width_die - 0.5*(full_refs.xmax - full_refs.xmin)
    full_refs.ymax = die.center[1] - distance_structure_from_center*width_die + 0.5*(full_refs.ymax - full_refs.ymin)        

    text_die = full_device << pg.text(text=name_of_chip, size=1.8e5, justify='left', layer=layer_text, font='DEPLOF')
    text_die.ymax = die.ymax - distance_text_from_border
    text_die.xmin = die.center[0] + distance_text_from_middle


    full_device = pg.union(full_device, by_layer = True);
    
    """
    Saving the design, along with it layer properties, and the parameters that 
    created it. We save the parameters using pickling.
    """

    if not os.path.isdir('Designs'):
        os.mkdir('Designs')
    
    pu.write_lyp(f'Designs/{name_of_chip}_layerinfo.lyp', layerset = lys)


    full_device.write_gds(filename = f'Designs/{name_of_chip}.gds', # Output GDS file name
                unit = 1e-9,                  # Base unit (1e-6 = microns)
                precision = 1e-9,             # Precision / resolution (1e-9 = nanometers)
                auto_rename = True,           # Automatically rename cells to avoid collisions
                max_cellname_length = 28,     # Max length of cell names
                cellname = 'toplevel'         # Name of output top-level cell
            )

def build_4_chip():
    
    print("|------------------| 1%")
    
    dot_params = inner_params['dot-params']
    amount_of_dots1 = dot_params['amount_of_dots_1']
    amount_of_dots2 = dot_params['amount_of_dots_2']
    amount_of_dots3 = dot_params['amount_of_dots_3']
    amount_of_dots4 = dot_params['amount_of_dots_4']
    
    print("|===---------------| 16.6%")
    
    local_marker_position_1 = inner_params['local_marker_positions_1']
    local_marker_position_2 = inner_params['local_marker_positions_2']
    local_marker_position_3 = inner_params['local_marker_positions_3']
    local_marker_position_4 = inner_params['local_marker_positions_4']

    print("|======------------| 33.3%")
    
    marker_structure_1, _, list_of_numbers_1, ext_ports_1, bond_numbers_1, layer_pads_1, prod_pads_1 = build_2xN(0, amount_of_dots1, layer_dot_array, layer_barrier_array, layer_screening_array, layer_ohmic_array, layers_outer, local_marker_position_1)
    marker_structure_2, _, list_of_numbers_2, ext_ports_2, bond_numbers_2, layer_pads_2, prod_pads_2 = build_2xN(1, amount_of_dots2, layer_dot_array, layer_barrier_array, layer_screening_array, layer_ohmic_array, layers_outer, local_marker_position_2)
    marker_structure_3, _, list_of_numbers_3, ext_ports_3, bond_numbers_3, layer_pads_3, prod_pads_3 = build_2xN(2, amount_of_dots3, layer_dot_array, layer_barrier_array, layer_screening_array, layer_ohmic_array, layers_outer, local_marker_position_3)
    marker_structure_4, _, list_of_numbers_4, ext_ports_4, bond_numbers_4, layer_pads_4, prod_pads_4 = build_2xN(3, amount_of_dots4, layer_dot_array, layer_barrier_array, layer_screening_array, layer_ohmic_array, layers_outer, local_marker_position_4)
    
    print("|=========---------| 50%")
    
    _, sem_structure_1, _, _, _, _, _ = build_2xN(0, amount_of_dots1, layer_dot_SEM_array, layer_barrier_SEM_array, layer_screening_SEM_array, layer_ohmic_SEM_array, layers_outer, local_marker_position_1)
    _, sem_structure_2, _, _, _, _, _ = build_2xN(1, amount_of_dots2, layer_dot_SEM_array, layer_barrier_SEM_array, layer_screening_SEM_array, layer_ohmic_SEM_array, layers_outer, local_marker_position_2)
    _, sem_structure_3, _, _, _, _, _ = build_2xN(2, amount_of_dots3, layer_dot_SEM_array, layer_barrier_SEM_array, layer_screening_SEM_array, layer_ohmic_SEM_array, layers_outer, local_marker_position_3)
    _, sem_structure_4, _, _, _, _, _ = build_2xN(3, amount_of_dots4, layer_dot_SEM_array, layer_barrier_SEM_array, layer_screening_SEM_array, layer_ohmic_SEM_array, layers_outer, local_marker_position_4)

    print("|============------| 66.6%")

    full_structure_1 = build_and_route_bond_pads(marker_structure_1, list_of_numbers_1, bond_numbers_1, ext_ports_1, prod_pads_1, layer_pads_1, width_bonding_box = width_bonding, nested = True)
    full_structure_2 = build_and_route_bond_pads(marker_structure_2, list_of_numbers_2, bond_numbers_2, ext_ports_2, prod_pads_2, layer_pads_2, width_bonding_box = width_bonding, nested = True)
    full_structure_3 = build_and_route_bond_pads(marker_structure_3, list_of_numbers_3, bond_numbers_3, ext_ports_3, prod_pads_3, layer_pads_3, width_bonding_box = width_bonding, nested = True)
    full_structure_4 = build_and_route_bond_pads(marker_structure_4, list_of_numbers_4, bond_numbers_4, ext_ports_4, prod_pads_4, layer_pads_4, width_bonding_box = width_bonding, nested = True)

    print("|===============---| 83.3%")

    structures_to_place = [full_structure_1, full_structure_2, full_structure_3, full_structure_4]

    sem_structures_to_place = [sem_structure_1, sem_structure_2, sem_structure_3, sem_structure_4]

    join_structures_on_chip(structures_to_place, sem_structures_to_place, layer_dummy_bonds)

    print("|==================| 100%")
    

def build_single_test():
    
    print("|------------------| 1%")
    
    dot_params = inner_params['dot-params']
    amount_of_dots1 = dot_params['amount_of_dots_1']
    
    print("|===---------------| 16.6%")
    
    local_marker_position_1 = inner_params['local_marker_positions_1']

    print("|======------------| 33.3%")
    
    marker_structure_1, _, list_of_numbers_1, ext_ports_1, bond_numbers_1, layer_pads_1, prod_pads_1 = build_2xN(0, amount_of_dots1, layer_dot_array, layer_barrier_array, layer_screening_array, layer_ohmic_array, layers_outer, local_marker_position_1)

    print("|=========---------| 50%")
    
    print("|============------| 66.6%")

    full_structure_1 = build_and_route_bond_pads(marker_structure_1, list_of_numbers_1, bond_numbers_1, ext_ports_1, prod_pads_1, layer_pads_1, width_bonding_box = width_bonding, nested = True)

    print("|===============---| 83.3%")

    structures_to_place = [full_structure_1]

    join_structures_test(structures_to_place)

    print("|==================| 100%")
    