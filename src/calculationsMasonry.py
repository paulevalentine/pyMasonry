import math

def pier_load(left_opening:float, pier_width:float, right_opening:float, wall_height:float, wall_load:float)->float:
    bearing_length = 0.150
    spread_angle = 60 * math.pi / 180

    # calculate the reactions from the lintels over the openings
    left_reaction = left_opening / 2 * wall_load
    right_reaction = right_opening / 2 * wall_load

    # calculate the lintel reaction spread lengths
    max_spread_length = bearing_length + wall_height / (2*math.tan(spread_angle))
    spread_length = min(max_spread_length, pier_width)

    # calculate the lintel UDLs
    left_udl = left_reaction / spread_length
    right_udl = right_reaction / spread_length

    # calculate the maximum UDL applied to the pier
    if 2 * spread_length > pier_width: # spread of lintel bearings overlaps
        pier_udl = (left_reaction + right_reaction) / pier_width + wall_load
    else:
        max_lintel_udl = max(left_reaction, right_reaction) / spread_length
        pier_udl = max_lintel_udl + wall_load

    print(f"The maximum UDL on the pier = {pier_udl:.2f}kN/m")
    return pier_udl


