import pandas as pd
import numpy as np
from geomeppy import IDF
from geomeppy import view_geometry
import pathlib
from matplotlib import pyplot as plt
import os
import subprocess # This should be used instead of os.system. Fix it later.
from geomeppy.geom.polygons import (
    break_polygons,
    Polygon2D,
    Polygon3D,
    Vector2D,
    Vector3D,
)
from geomeppy.geom.transformations import align_face, invert_align_face
"""
Use
1) Make a poly from an IDF surface with the "IDFsurf_to_poly" func
2) Trans the poly to x-y plane like this: "poly_trans = align_face(poly).order_points('upperleftcorner')" and then make it to 2d with the "project_to_2D" function
3) Use the (still 3d) but already translated poly to make a bbox and project it to 2d: "poly_trans_bbox_2d = poly_trans.bounding_box.project_to_2D()"
4) Make a test grid in 2D using the bbox we made and the "grid_2d" func
5) Make a list of test points in the poly in 2d with the "gridpoints_in_poly_2d" func
6) Use the original poly, the surface name the list of testpoints we made with the "create_sensor_points" func to make a dict, that has the sensorpoint name, the sensorpoints for that IDF suface, the direction vectors of the sensor points
7) Make a list of the previous stuff made with (6) and translate it to daysim .pts file and a "keep track of sensorpoints" .kts supplementary file with the "translate_to_ds_pts" func
"""
def IDFsurf_to_poly(surface):
    "Makes a poly from IDF surface"
    # get its vertices
    corstrings = ['_Xcoordinate', '_Ycoordinate', '_Zcoordinate']
    vertnums = list(range(1, int(surface['Number_of_Vertices'])+1, 1))
    surf_verts = []
    for vn in vertnums:
        vert_cords = tuple([surface['Vertex_'+ str(vn) + cs] for cs in corstrings])
        surf_verts.append(vert_cords)
    # make a 3d polygon from vertices
    poly = Polygon3D(surf_verts)
    return poly

def grid_2d_depr1(bbox_2d, n):
    """Creates a grid over the bbox with sides divided to n. TODO: n-ratio"""
    xs = (bbox_2d.vertices_list[0][0], bbox_2d.vertices_list[2][0])
    ys = (bbox_2d.vertices_list[0][1], bbox_2d.vertices_list[2][1])
    xv = np.linspace(xs[0], xs[1], n)
    yv = np.linspace(ys[0], ys[1], n) # this should be a func of the reatio of width height of the bbox
    grid = np.meshgrid(xv, yv)
    return grid

def grid_2d(bbox_2d, d):
    """Creates a grid over the bbox with sided divided to n. TODO: n-ratio"""
    a = abs((bbox_2d[0] - bbox_2d[2])[0])
    b = abs((bbox_2d[0] - bbox_2d[2])[1])
    n = d * a * b
    #y = (n/(b/a + 1)) - 1
    #x = (b/a) * y
    y1 = (-(a/b + 1) + np.sqrt((a/b + 1)**2 - 4*a/b*(1 - a*b*d))) / (2 * a/b)
    y2 = (-(a/b + 1) - np.sqrt((a/b + 1)**2 - 4*a/b*(1 - a*b*d))) / (2 * a/b)
    x1 = y1* a/b
    x2 = y2* a/b
    y = max(y1, y2)
    x = y * a/b
    print('a={}; b={}; x={}; y={}; n={}'.format(a,b,x,y,n))
    xs = (bbox_2d.vertices_list[0][0], bbox_2d.vertices_list[2][0])
    ys = (bbox_2d.vertices_list[0][1], bbox_2d.vertices_list[2][1])
    xv = np.linspace(xs[0], xs[1], round(x+1, 0))
    yv = np.linspace(ys[0], ys[1], round(y+1, 0)) # this should be a func of the reatio of width height of the bbox
    #xv = np.linspace(xs[0], xs[1], a/x, 0)
    #yv = np.linspace(ys[0], ys[1], b/y, 0) # this should be a func of the reatio of width height of the bbox
    grid = np.meshgrid(xv, yv)
    return grid

def grid_2d_pvmodule(n_row, n_col, bbox_2d):
    "This makes a rectangular grid over the pv cells. Use the gridpoints_in_poly_2d on the list of arrays it returns, because this produces some extra points as well, that are needed to be clipped off"
    d = (n_col*n_row) + (n_row + 1) + (n_col + 1) - 1 # calc the point density to have 1 sp over each cell
    test_grid = grid_2d(bbox_2d=bbox_2d, d=d) # make a 2D test grid on the x-y plane (daypym) This grid contains the edges as well. Not good yet. We will push this with half cell-size in eac direction.
    lr = test_grid[0].max() - test_grid[0].min() # calc the length of the module in the row-direction
    lc = test_grid[1].max() - test_grid[1].min() # calc the length of the module in the column-direction
    dr = (lr/n_row) / 2
    dc = (lc/n_col) / 2
    test_grid_transp = [test_grid[0] + dr, test_grid[1] + dc] # Pushing the grid to the center of the cells. Now we have a row and column of points off the surface of the pv. These are needed to be removed with gridpoints_in_poly_2d func!
    return test_grid_transp

def point_in_poly_depr1(point, poly):
    "True if a 2D point is in a 2D poly. Deprecated: This left out the points on the right hand side edges" # add dynamic quasi infinite point distance with bounding box clue
    qip1 = Vector2D(10000, 0) # quasi-infinite point
    qip2 = Vector2D(10000, 1) # quasi-infinite point 2
    test_poly = Polygon2D([point, qip1, qip2])#.order_points('upperleftcorner') # maybe we wouldn't need this anyway
    #poly = poly.order_points('upperleftcorner') # maybe we wouldn't need this anyway
    inter = poly.intersect(test_poly)
    if len(inter) == 0:
        return False
    else:
        no_inter = len(inter[0])
        #print(inter)
        return True if no_inter == 3 else False

def point_in_poly(point, poly):
    "True if a 2D point is in or on the edge of a 2D poly."
    # TODO: add dynamic quasi infinite point distance with bounding box clue
    qip1_r = Vector2D(10000, 0) # quasi-infinite point 1 to the right
    qip2_r = Vector2D(10000, 1) # quasi-infinite point 2 to the right
    qip1_l = Vector2D(-10000, 0) # quasi-infinite point 1 to the left
    qip2_l = Vector2D(-10000, 1) # quasi-infinite point 2 to the left
    test_poly_r = Polygon2D([point, qip1_r, qip2_r])#.order_points('upperleftcorner') # maybe we wouldn't need this anyway
    test_poly_l = Polygon2D([point, qip1_l, qip2_l])#.order_points('upperleftcorner') # maybe we wouldn't need this anyway
    #poly = poly.order_points('upperleftcorner') # maybe we wouldn't need this anyway
    inter_r = poly.intersect(test_poly_r)
    inter_l = poly.intersect(test_poly_l)
    if len(inter_r) == 0:
        inter_r.append([0])
    if len(inter_l) == 0:
        inter_l.append([0])
    no_inter_r = len(inter_r[0])
    no_inter_l = len(inter_l[0])
    return True if no_inter_r == 3 or no_inter_l == 3 else False # we will see how robust this is # check this again later

def gridpoints_in_poly_2d(grid_2d, poly_2d):
    "Returns a list of points, that are inside a 2d poly. Takes a meshgrid (made with grid_2d) and a 2d poly as input"
    pip = [] # list of points in poly
    for x, y in list(zip(grid_2d[0], grid_2d[1])):
        for cx, cy in list(zip(x, y)):
            test_point = Vector2D(cx, cy)
            if point_in_poly(point=test_point, poly=poly_2d):
                pip.append(test_point)
    return pip

def gridpoints_in_poly_2d2(grid_2d, poly_2d):
    "Returns a list of points, that are inside a 2d poly. Takes a meshgrid (made with grid_2d) and a 2d poly as input"
    pip = [] # list of points in poly
    for x, y in list(zip(grid_2d[0], grid_2d[1])):
        for cx, cy in list(zip(x, y)):
            test_point = Vector2D(cx, cy)
            if point_in_poly2(point=test_point, poly=poly_2d):
                pip.append(test_point)
    return pip

def create_sensor_points(surf_name, points_in_poly_2d, original_poly, sp_offset=0.01, sp_pos_round=3):
    "Returns a dict {'surf_name':surface name, 'sensor_points':list of sensor points, 'sp_ori':list of normal vectors}. sesor points are translated with sp_offset from the original plane. Takes points_in_poly_2d (made with gridpoints_in_poly_2d function) and the original poly as input. We can round the coords of the sps. Make sure, that the rounding has more digits, than the smallest distances and the offset. Default is 1 cm offset and 1 mm rounding."
    # make a 3D polygon from these points and rotate them back
    sensor_points_poly_trans = Polygon3D(points_in_poly_2d)
    sensor_points_poly = invert_align_face(original=original_poly, poly2=sensor_points_poly_trans)
    # give the sensor point poly a nudge outwards (let's say 1 cm abowe the surface)
    #####just add the 1/100th of the normal vector to the points
    #normal_vectors = original_poly.normal_vector.as_array() # for some reason as_array returns a vector, that has a good direction, but the length is not 1
    normal_vectors = np.array([original_poly.normal_vector[0], original_poly.normal_vector[1], original_poly.normal_vector[2]])
    transl_vect = normal_vectors * sp_offset
    sensor_points = [list(np.add(tuple([round(sp[i], sp_pos_round) for i in range(len(sp))]), transl_vect)) for sp in sensor_points_poly.vertices_list] # sorry
    sensor_points_surf = {'surf_name':surf_name, 'surf_coords':original_poly.vertices_list, 'sensor_points':sensor_points, 'sp_ori':list(normal_vectors)} # make a dict for the sensor points with name points and vectors, later export this to DS file
    return sensor_points_surf

def translate_to_ds_pts(surf_sensor_points, outputpath_pts, outputpath_kts):
    x, y, z, vx, vy, vz, surf_name = ([] for i in list(range(7)))
    for surf in surf_sensor_points:
        for i in list(range(len(surf['sensor_points']))):
            x.append(surf['sensor_points'][i][0])
            y.append(surf['sensor_points'][i][1])
            z.append(surf['sensor_points'][i][2])
            vx.append(surf['sp_ori'][0])
            vy.append(surf['sp_ori'][1])
            vz.append(surf['sp_ori'][2])
            surf_name.append(surf['surf_name'])
    spdf = pd.DataFrame(data={'x':x, 'y':y, 'z':z, 'vx':vx, 'vy':vy, 'vz':vz, 'surf_name':surf_name})
    ptsfile = open(outputpath_pts, 'w')
    ktsfile = open(outputpath_kts, 'w')
    for p in spdf.index:
        ptsfile.write(str(spdf['x'][p]) +' '+ str(spdf['y'][p]) +' '+ str(spdf['z'][p]) +' '+ str(spdf['vx'][p]) +' '+ str(spdf['vy'][p]) +' '+ str(spdf['vz'][p]) + '\n')
        ktsfile.write(str(spdf['surf_name'][p] + '\n'))
    ptsfile.close()
    ktsfile.close()
    #print('translate_to_ds_pts: .pts and .kts files saved')
    return

def view_idf_to_ax(fname=None, idf_txt=None, test=False):
    """This is originally from https://github.com/jamiebull1/geomeppy/blob/master/geomeppy/view_geometry.py
    This just returns an ax instead of viewing it on order to  plot it together with the sensorpoints"""
    from geomeppy.view_geometry import _get_collection, _get_collections, _get_surfaces, _get_limits # probably these should not be imported here
    from io import StringIO
    from eppy.iddcurrent import iddcurrent
    # type: (Optional[str], Optional[str], Optional[bool]) -> None

    if fname and idf_txt:
        raise ValueError("Pass either fname or idf_txt, not both.")
    # set the IDD for the version of EnergyPlus
    iddfhandle = StringIO(iddcurrent.iddtxt)
    if IDF.getiddname() is None:
        IDF.setiddname(iddfhandle)

    if fname:
        # import the IDF
        idf = IDF(fname)
    elif idf_txt:
        idf = IDF()
        idf.initreadtxt(idf_txt)
    # create the figure and add the surfaces
    ax = plt.axes(projection="3d")
    collections = _get_collections(idf, opacity=0.5)
    for c in collections:
        ax.add_collection3d(c)

    # calculate and set the axis limits
    limits = _get_limits(idf=idf)
    ax.set_xlim(limits["x"])
    ax.set_ylim(limits["y"])
    ax.set_zlim(limits["z"])
    return ax

def view_idf_and_sps(p_name, idf_name, sps):
    "To view the e+ IDF and the DS sensorpoints together. Utilizing modified version of Geomeppy view_idf function: https://github.com/jamiebull1/geomeppy/"
    # TODO implement auto save fig: Save_fig False: no saving, True: save it to /geo
    surfcoords = []
    polys = []
    for surf in sps:
        polys.append(Polygon3D(surf['surf_coords']))
        for sp in surf['sensor_points']:
            surfcoords.append((sp[0], sp[1], sp[2]))
    xs = [c[0] for c in surfcoords]
    ys = [c[1] for c in surfcoords]
    zs = [c[2] for c in surfcoords]
    ax2 = view_idf_to_ax(fname=idf_name, idf_txt=None, test=False)
    ax2.scatter(xs, ys, zs, marker='o', s=2, c='k')
    try:
        plt.savefig('geo/{}.png'.format(p_name))
    except:
        plt.savefig('{}.png'.format(p_name))
    plt.show(block=False)

def ds_epw2wea(input_epw, output_wea, auto_path=True, return_wea_info=True):
    """execute DS program: epw2wea <input_file_name.epw> <output_file_name.wea>. If autopath not True, use full paths. With autopath it will use the pre-defined folder structure"""
    if auto_path == True:
        cstr = r'epw2wea {} wea\{}.wea'.format(input_epw, output_wea)
    else:
        cstr = r'epw2wea {} {}'.format(input_epw, output_wea)
    os.system(cstr)
    if return_wea_info == True:
        #result = subprocess.run([cstr], stdout=subprocess.PIPE)
        result = subprocess.check_output(cstr, stderr=subprocess.STDOUT)
        #print(stderr)
        result = result.decode('utf-8').split('\r\n')
        kl = [i.split()[0] for i in result[:-1]]
        vl = [i.split()[1] for i in result[:-1]]
        return dict(list(zip(kl, vl)))

def rad_obj2rad(input_obj, output_rad, auto_path=True):
    """execute DS program: obj2rad exported_geom.obj > test3_geom.rad. If autopath not True, use full paths. With autopath it will use the pre-defined folder structure"""
    if auto_path == True:
        cstr = r'obj2rad geo\{} > geo\{}.rad'.format(input_obj, output_rad)
    else:
        cstr = r'obj2rad {} > {}'.format(input_obj, output_rad)
    os.system(cstr)

def write_ds_hea(p_name, p_dir, bin_dir, uni_mat_file, site_info, model_info=None, radiance_params=None):
    """Write the hea file to the project dir"""
    if radiance_params == None:
        radiance_params = {  # check later these default params, to set to something practical for irrad sim
                           'ab':2,
                           'ad':1500,
                           'as':20,
                           'ar':300,
                           'aa':0.05,
                           'lr':6,
                           'st':0.15,
                           'sj':1,
                           'lw':0.004,
                           'dj':0,
                           'ds':0.2,
                           'dr':2,
                           'dp':512
                          }
    if model_info == None:
        # get num of sensorpoints
        with open(p_dir + '\pts\\' + p_name + '.pts') as f:
            for n, l in enumerate(f):
                pass
        n = n + 1
        model_info = {
                      'radiance_source_files':'2, {}\\geo\{}.rad, {}'.format(p_dir, p_name, uni_mat_file),
                      'material_file':'\\geo\ds_mat.rad', # this is output made by daysim
                      'geometry_file':'\\geo\ds_geo.rad', # this is output made by daysim
                      'sensor_file':'\\pts\\' + p_name + '.pts', # this is input made by us before
                      'shading':'1 {} \\res\{}.dc \\res\{}.ill'.format(p_name, p_name, p_name), # this is output made by daysim
                      'sensor_file_unit':'2 ' * n,
                      #'DC_file_format':'Daysim_original',
                      #'DC_file_format':'dds',
                      'DDS_sensor_file':r'\\res\\{}.sen'.format(p_name),
	                  'DDS_file':r'\\res\\{}.dds'.format(p_name)
                     }

    project_info = {
                    'project_name':p_name,
                    'project_directory':p_dir,
                    'bin_directory':bin_dir,
                    'temp_directory':p_dir + r'/tmp/',
                    #'material_directory':mat_dir
                    }
    # writing the hea file
    hea = open(p_dir +r'\\'+ p_name + r'.hea', 'w')
    hea.write('# DAYSIM 4 project file created with daypym' + '\n\n') # 'title'
    for key, value in project_info.items():
        hea.write(key + '\t\t' + str(value) + '\n')
    hea.write('\n')
    hea.write('#### SITE INFORMATION ####\n')
    for key, value in site_info.items():
        hea.write(key + '\t\t' + str(value) + '\n')
    hea.write('\n')
    hea.write('#### MODEL INFORMATION ####\n')
    for key, value in model_info.items():
        hea.write(key + '\t\t' + str(value) + '\n')
    hea.write('\n')
    hea.write('#### RADIANCE PARAMETERS ####\n')
    for key, value in radiance_params.items():
        hea.write(key + '\t\t' + str(value) + '\n')
    hea.close()

def ds_radfiles2daysim(p_name, hea_file_name=None):
    """Run radfiles2daysim. Use project name for automated workflow or the hea file name (with extension) for manual workflow."""
    if hea_file_name == None:
        cstr = r'radfiles2daysim {}.hea -m -g'.format(p_name)
    else:
        cstr = r'radfiles2daysim {} -m -g'.format(hea_file_name)
    os.system(cstr)

def ds_ds_shortterm(p_name, hea_file_name=None):
    """Run gen_dc. Use project name for automated workflow or the hea file name (with extension) for manual workflow. New short timestep is in the .hea file, given via the site_info."""
    if hea_file_name == None:
        cstr = r'ds_shortterm {}.hea'.format(p_name)
    else:
        cstr = r'ds_shortterm {}.hea'.format(hea_file_name)
    os.system(cstr)

def ds_gen_dc(p_name, hea_file_name=None):
    """Run gen_dc. Use project name for automated workflow or the hea file name (with extension) for manual workflow."""
    if hea_file_name == None:
        cstr = r'gen_dc {}.hea -dds'.format(p_name)
    else:
        cstr = r'gen_dc {} -dds'.format(hea_file_name)
    os.system(cstr)

def ds_ds_illum(p_name, hea_file_name=None):
    """Run ds_illum. Use project name for automated workflow or the hea file name (with extension) for manual workflow."""
    if hea_file_name == None:
        cstr = r'ds_illum {}.hea -dds'.format(p_name)
    else:
        cstr = r'ds_illum {} -dds'.format(hea_file_name)
    os.system(cstr)
