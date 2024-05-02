#Function backup, the real one is in the file 'export_topo_operative.py'

def export_flac3D_input(geo_model, path=None, filename='geomodel.f3grid'):
 
    # create vertices and elements
    vertices, elements, groups = __build_vertices_elements_groups__(geo_model)

    # open output file
    # if not path:
    #    path = './'
    # if not os.path.exists(path):
    #    os.makedirs(path)

    out = open(path + filename, 'w')

    # write gridpoints
    out.write("*GRIDPOINTS")
    for i, vertice in enumerate(vertices):
        out.write(f"\nG {i + 1} {vertice[0]} {vertice[1]} {vertice[2]}")

    # write elements
    out.write('\n*ZONES')
    zone_counter = 1 

    new_column_order = [1, 2, 4, 5, 3, 8, 6, 7]

    for i, elem in enumerate(elements):
        out.write(f'\nZ B8 {zone_counter}')  
        for x in new_column_order:
            out.write(f" {elem[x]}")
        zone_counter += 1 

    # make groups
    out.write('\n*GROUPS\n')
    for grp_name, grp in groups.items():
        out.write(f'ZGROUP \"{grp_name}\"\n')
        count = 0
        for x in grp:
            out.write(f"{x} ")
            count += 1
            if count == 8:
                out.write("\n")
                count = 0
        if count != 0: out.write("\n")

    out.close()
    print("Successfully exported geological model as FLAC3D input to " + path)
    return 


def __build_vertices_elements_groups__(geo_model):

    # get model information
    nx, ny, nz = geo_model.grid.regular_grid.resolution
    xmin, xmax, ymin, ymax, zmin, zmax = geo_model.solutions.grid.regular_grid.extent

    # create vertices array
    dx, dy, dz = (xmax - xmin) / nx, (ymax - ymin) / ny, (zmax - zmin) / nz
    n_vertices = (nx + 1) * (ny + 1) * (nz + 1)
    vertices = np.zeros((n_vertices, 3), dtype='f8')
    vertices_ids = np.arange(n_vertices)  # used to generate coordinate
    vertices[:, 0] = vertices_ids % (nx + 1) * dx + xmin
    vertices[:, 1] = (vertices_ids % ((nx + 1) * (ny + 1))) // (nx + 1) * dy + ymin
    vertices[:, 2] = vertices_ids // ((nx + 1) * (ny + 1)) * dz + zmin

    # build elements
    n_elements = nx * ny * nz
    element_ids = np.arange(n_elements)  # used to generate elems
    elements = np.zeros((n_elements, 9), dtype='i8')
    i = element_ids % nz
    j = element_ids // nz % ny
    k = element_ids // (nz * ny)
    elements[:, 0] = 8  # all hex
    elements[:, 1] = 1 + i * (nx + 1) * (ny + 1) + j * (nx + 1) + k
    elements[:, 2] = elements[:, 1] + 1
    elements[:, 3] = elements[:, 2] + (nx + 1)
    elements[:, 4] = elements[:, 3] - 1
    elements[:, 5] = elements[:, 1] + ((nx + 1) * (ny + 1))
    elements[:, 6] = elements[:, 5] + 1
    elements[:, 7] = elements[:, 6] + (nx + 1)
    elements[:, 8] = elements[:, 7] - 1

    # build groups
    lith_ids = np.round(geo_model.solutions.lith_block)
    lith_ids = lith_ids.astype(int)
    sids = dict(zip(geo_model._surfaces.df['surface'], geo_model._surfaces.df['id']))
    groups = {}
    for region_name, region_id in sids.items():
        cell_ids = np.where(lith_ids == region_id)[0] + 1
        if not len(cell_ids): continue
        groups[region_name] = cell_ids

    # remove element above topography
    #mask_topo = geo_model._grid.regular_grid.mask_topo
    shape = geo_model._grid.regular_grid.resolution
    shape_tot = shape[0]*shape[1]*shape[2]
    mask_topo = geo_model._grid.regular_grid.mask_topo
    if(mask_topo.size > 0):
        inactive_cells = mask_topo.reshape(shape_tot)
    else:
        inactive_cells = None
    if np.any(inactive_cells):
        # update correspondance
        new_id_vertices = np.zeros(len(vertices), dtype='i8')
        new_id_elements = np.zeros(len(elements), dtype='i8')
        # remove inactive cell
        elements = elements[~inactive_cells]
        new_id_elements[~inactive_cells] = np.arange(len(elements)) + 1
        # remove deleted vertices
        cond = np.isin(np.arange(len(vertices)) + 1, elements[:, 1:].flatten())
        vertices = vertices[cond]
        new_id_vertices[cond] = np.arange(len(vertices)) + 1
        # renumber vertices in element
        elements[1:] = new_id_vertices[elements[1:] - 1]
        # renumber groups
        for grp in groups.values():
            grp[:] = new_id_elements[grp - 1]

    return vertices, elements, groups