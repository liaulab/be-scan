from pymol import cmd, stored, math
	
def loadBfacts(
    mol, startaa=1, source="newBfactors.txt", 
    visual="white_red", legend_name="count",
    scale=1, lower=None, upper=None, radius=0.2,
    color_names=["myred", "myblue"],
    color_codes=["#6666FF", "#FF6666"],
    ):

    """
    Replaces B-factors with a list of values contained in a plain txt file. 
    The plain text file should contain only a list of values of the same length as the protein of interest. 
    """

    # COLOR SCALE #
    for color_name, color_code in zip(color_names, color_codes):
        r = int(color_code[1:3], 16) / 255.0
        g = int(color_code[3:5], 16) / 255.0
        b = int(color_code[5:7], 16) / 255.0
        cmd.set_color(color_name, [r, g, b])

    # PLOT DATA #
    scale = float(scale)
    obj=cmd.get_object_list(mol)[0]
    cmd.alter(mol,"b=-1.0")
    inFile = open(source, 'r')
    counter=int(startaa)
    bfacts=[]
    for line in inFile.readlines():
        bfact=float(line) * scale
        bfacts.append(bfact)
        cmd.alter("%s and resi %s and n. CA"%(mol,counter), "b=%s"%bfact)
        counter=counter+1

    if lower is not None: bf_min = float(lower)
    else: bf_min = min(bfacts)
    if upper is not None: bf_max = float(upper)
    else: bf_max = max(bfacts)

    # COLOR CARTOON #
    cmd.show_as("cartoon", mol)
    cmd.cartoon("putty", mol)
    cmd.set("cartoon_putty_scale_min", bf_min, obj)
    cmd.set("cartoon_putty_scale_max", bf_max, obj)
    cmd.set("cartoon_putty_transform", 0, obj)
    cmd.set("cartoon_putty_radius", radius, obj)

    # CHOOSE COLOR PALETTE #
    if visual == "rainbow":
        cmd.spectrum("b", "rainbow", "%s and n. CA " %mol, minimum=bf_min, maximum=bf_max)
    #     cmd.ramp_new(legend_name, obj, [scale*bf_min, scale*bf_max], "rainbow")

    # if visual=="blue_red":
    #     cmd.spectrum("b", "myblue myred", "%s and n. CA " %mol, minimum=bf_min, maximum=bf_max)
    #     cmd.ramp_new(legend_name, obj, [0, scale*bf_max], ["myblue", "myred"])

    if visual=="white_red":
        print(bf_min, bf_max)
        cmd.spectrum("b", "white red", "%s and n. CA " %mol, minimum=bf_min, maximum=bf_max)
    #     cmd.ramp_new(legend_name, obj, [0, bf_max], ["white", "red"])

    cmd.recolor()

cmd.extend("loadBfacts", loadBfacts);
