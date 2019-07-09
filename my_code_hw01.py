#-- my_code_hw01.py
#-- hw01 GEO1015/2018
#-- Yifang Zhao
#-- 4798899
#-- Jinglan Li
#-- 4781937


import math
import scipy.spatial
import numpy


def nn_interpolation(list_pts_3d, j_nn):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_nn:        the parameters of the input for "nn"
    Output:
        returns the value of the area
 
    """  
    print("=== Nearest neighbour interpolation ===")

    # print("cellsize:", j_nn['cellsize'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # d, i = kd.query(p, k=1)
    
    cellsize = j_nn['cellsize'];
    xcoord = []
    ycoord = []
    list_pts_2d = []
    for p in list_pts_3d:
        xcoord.append(p[0])
        ycoord.append(p[1])
        list_pts_2d.append([p[0], p[1]])
    xmin = min(xcoord)
    xmax = max(xcoord)
    ymin = min(ycoord)
    ymax = max(ycoord)
    x = xmin
    y = ymin
    rows = math.ceil((ymax - ymin) / cellsize)
    columns = math.ceil((xmax - xmin) / cellsize)
    
    kd = scipy.spatial.KDTree(list_pts_2d)
    
    hull = scipy.spatial.ConvexHull(list_pts_2d)
    list_vertices = []
    for i in hull.vertices:
        list_vertices.append(list_pts_2d[i])
    list_vertices.append(list_vertices[0])
    
    with open(j_nn['output-file'], 'w') as f:
        f.write('NCOLS ' + str(columns) + '\n')
        f.write('NROWS ' + str(rows) + '\n')
        f.write('XLLCORNER ' + str(xmin) + '\n')
        f.write('YLLCORNER ' + str(ymin) + '\n')
        f.write('CELLSIZE ' + str(cellsize) + '\n')
        f.write('NODATA_VALUE -9999' + '\n')
        for row in range(rows, 0, -1):  #from the left-upper corner
            for column in range(columns):
                p = [x + cellsize * column + cellsize / 2, y + cellsize * row - cellsize / 2]
                if is_point_inside_single_polygon(p, list_vertices):
                    d, i = kd.query(p, k=1)
                    p_nn = list_pts_3d[i]
                    f.write(str(p_nn[2]))
                else:
                    f.write('-9999')
                if column != columns - 1:
                    f.write(" ")
            f.write("\n")
    print("File written to", j_nn['output-file'])



def idw_interpolation(list_pts_3d, j_idw):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with IDW
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_idw:       the parameters of the input for "idw"
    Output:
        returns the value of the area
 
    """  
    print("=== IDW interpolation ===")

    # print("cellsize:", j_idw['cellsize'])
    # print("radius:", j_idw['radius'])

    #-- to speed up the nearest neighbour us a kd-tree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
    # kd = scipy.spatial.KDTree(list_pts)
    # i = kd.query_ball_point(p, radius)
    
    cellsize = j_idw['cellsize'];
    radius = j_idw['radius']
    xcoord = []
    ycoord = []
    list_pts_2d = []
    for p in list_pts_3d:
        xcoord.append(p[0])
        ycoord.append(p[1])
        list_pts_2d.append([p[0], p[1]])
    xmin = min(xcoord)
    xmax = max(xcoord)
    ymin = min(ycoord)
    ymax = max(ycoord)
    x = xmin
    y = ymin
    rows = math.ceil((ymax - ymin) / cellsize)
    columns = math.ceil((xmax - xmin) / cellsize)
    
    kd = scipy.spatial.KDTree(list_pts_2d)
    
    hull = scipy.spatial.ConvexHull(list_pts_2d)
    list_vertices = []
    for i in hull.vertices:
        list_vertices.append(list_pts_2d[i])
    list_vertices.append(list_vertices[0])
    
    with open(j_idw['output-file'], 'w') as f:
        f.write('NCOLS ' + str(columns) + '\n')
        f.write('NROWS ' + str(rows) + '\n')
        f.write('XLLCORNER ' + str(xmin) + '\n')
        f.write('YLLCORNER ' + str(ymin) + '\n')
        f.write('CELLSIZE ' + str(cellsize) + '\n')
        f.write('NODATA_VALUE -9999' + '\n')
        for row in range(rows, 0, -1):  #from the left-upper corner
            for column in range(columns):
                p = [x + cellsize * column + cellsize / 2, y + cellsize * row - cellsize / 2]
                if is_point_inside_single_polygon(p, list_vertices):
                    height = 0
                    array_d, array_i = kd.query(p, k = 1000, distance_upper_bound = radius)
                    cutoff_idx = list(array_i).index(len(list_pts_2d))
                    if cutoff_idx == 0: # if no points were found inside the search area, use nn then.
                        d, i = kd.query(p, k=1)
                        height = list_pts_3d[i][2]
                    else:
                        num = 0
                        for j in range(cutoff_idx):
                            num += 1 / array_d[j]**2
                        for j in range(cutoff_idx):
                            height += (1/array_d[j]**2) / num * list_pts_3d[array_i[j]][2]
                    f.write(str(height))
                else:
                    f.write('-9999')
                if column != columns - 1:
                    f.write(" ")
            f.write("\n")
    
    print("File written to", j_idw['output-file'])


def tin_interpolation(list_pts_3d, j_tin):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with linear in TIN interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_tin:       the parameters of the input for "tin"
    Output:
        returns the value of the area
 
    """  
    print("=== TIN interpolation ===")

    #-- example to construct the DT
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html#scipy.spatial.Delaunay
    # dt = scipy.spatial.Delaunay([])
    
    cellsize = j_tin['cellsize'];
    xcoord = []
    ycoord = []
    list_pts_2d = []
    for p in list_pts_3d:
        xcoord.append(p[0])
        ycoord.append(p[1])
        list_pts_2d.append([p[0], p[1]])
    xmin = min(xcoord)
    xmax = max(xcoord)
    ymin = min(ycoord)
    ymax = max(ycoord)
    x = xmin
    y = ymin
    rows = math.ceil((ymax - ymin) / cellsize)
    columns = math.ceil((xmax - xmin) / cellsize)
    
    dt = scipy.spatial.Delaunay(list_pts_2d)
    
    hull = scipy.spatial.ConvexHull(list_pts_2d)
    list_vertices = []
    for i in hull.vertices:
        list_vertices.append(list_pts_2d[i])
    list_vertices.append(list_vertices[0])
    
    with open(j_tin['output-file'], 'w') as f:
        f.write('NCOLS ' + str(columns) + '\n')
        f.write('NROWS ' + str(rows) + '\n')
        f.write('XLLCORNER ' + str(xmin) + '\n')
        f.write('YLLCORNER ' + str(ymin) + '\n')
        f.write('CELLSIZE ' + str(cellsize) + '\n')
        f.write('NODATA_VALUE -9999' + '\n')
        for row in range(rows, 0, -1):  #from the left-upper corner
            for column in range(columns):
                p = [x + cellsize * column + cellsize / 2, y + cellsize * row - cellsize / 2]
                if is_point_inside_single_polygon(p, list_vertices):
##                    A = numpy.mat(numpy.array(list_pts_3d)[dt.simplices[dt.find_simplex(p)]])
##                    co = A.I * numpy.mat([1, 1, 1]).T
##                    height = (1 - co[0, 0]*p[0] - co[1, 0]*p[1]) / co[2, 0]
                    vtxs = numpy.array(list_pts_3d)[dt.simplices[dt.find_simplex(p)]]
                    vx = vtxs[:,0]
                    vy = vtxs[:,1]
                    vz = vtxs[:,2]
                    w1 = abs((1/2)*(p[0]*vy[1]+vx[1]*vy[2]+vx[2]*p[1]-p[0]*vy[2]-vx[1]*p[1]-vx[2]*vy[1]))
                    w2 = abs((1/2)*(p[0]*vy[2]+vx[2]*vy[0]+vx[0]*p[1]-p[0]*vy[0]-vx[2]*p[1]-vx[0]*vy[2]))
                    w3 = abs((1/2)*(p[0]*vy[0]+vx[0]*vy[1]+vx[1]*p[1]-p[0]*vy[1]-vx[0]*p[1]-vx[1]*vy[0]))
                    height = (vz[0]*w1 + vz[1]*w2 + vz[2]*w3)/(w1 + w2 + w3)
                    f.write(str(height))
                else:
                    f.write('-9999')
                if column != columns - 1:
                    f.write(" ")
            f.write("\n")
    
    print("File written to", j_tin['output-file'])


def kriging_interpolation(list_pts_3d, j_kriging):
    """
    !!! TO BE COMPLETED !!!
     
    Function that writes the output raster with ordinary kriging interpolation
     
    Input:
        list_pts_3d: the list of the input points (in 3D)
        j_kriging:       the parameters of the input for "kriging"
    Output:
        returns the value of the area
 
    """  
    print("=== Ordinary kriging interpolation ===")

    # Remove duplicate points
    clean_points_list = []
    for point1 in list_pts_3d:
    	repeated = False
    	for point2 in clean_points_list:
    		if point1[0] == point2[0] and point1[1] == point2[1]:
    			repeated = True
    	if repeated == False:
    		clean_points_list.append(point1)
    	else:
    		print("Repeated point: " + str(point1[0]) + " " + str(point1[1]))
    points_3d = numpy.array(clean_points_list)
    
    cellsize = j_kriging['cellsize'];
    radius = j_kriging['radius']
    xcoord = []
    ycoord = []
    points_2d = []
    for p in points_3d:
        xcoord.append(p[0])
        ycoord.append(p[1])
        points_2d.append([p[0], p[1]])
    xmin = min(xcoord)
    xmax = max(xcoord)
    ymin = min(ycoord)
    ymax = max(ycoord)
    x = xmin
    y = ymin
    rows = math.ceil((ymax - ymin) / cellsize)
    columns = math.ceil((xmax - xmin) / cellsize)
    
    kd = scipy.spatial.KDTree(points_2d)
    
    hull = scipy.spatial.ConvexHull(points_2d)
    list_vertices = []
    for i in hull.vertices:
        list_vertices.append(points_2d[i])
    list_vertices.append(list_vertices[0])
    
    with open(j_kriging['output-file'], 'w') as f:
        f.write('NCOLS ' + str(columns) + '\n')
        f.write('NROWS ' + str(rows) + '\n')
        f.write('XLLCORNER ' + str(xmin) + '\n')
        f.write('YLLCORNER ' + str(ymin) + '\n')
        f.write('CELLSIZE ' + str(cellsize) + '\n')
        f.write('NODATA_VALUE -9999' + '\n')
        for row in range(rows, 0, -1):  #from the left-upper corner
            for column in range(columns):
                p = [x + cellsize * column + cellsize / 2, y + cellsize * row - cellsize / 2]
                if is_point_inside_single_polygon(p, list_vertices):
                    array_d, array_i = kd.query(p, k = 1000, distance_upper_bound = radius)
                    cutoff_idx = list(array_i).index(len(points_2d))
                    if cutoff_idx == 0: # if no points were found inside the search area, use nn then.
                        d, i = kd.query(p, k=1)
                        height = list_pts_3d[i][2]
                    else:
                        A = numpy.mat(numpy.ones((cutoff_idx + 1, cutoff_idx + 1)))
                        d = numpy.mat(numpy.ones((cutoff_idx + 1, 1)))
                        for m in range(cutoff_idx):
                            d[m, 0] = g(dis(p, points_2d[array_i[m]]))
                            for n in range(cutoff_idx):
                                A[m, n] = g(dis(points_2d[array_i[m]], points_2d[array_i[n]]))
                        A[cutoff_idx, cutoff_idx] = 0
                        w = A.I * d
                        height = 0
                        for j in range(cutoff_idx):
                            height += w[j, 0] * points_3d[array_i[j]][2]
                    f.write(str(height))
                else:
                    f.write('-9999')
                if column != columns - 1:
                    f.write(" ")
            f.write("\n")
    
    print("File written to", j_kriging['output-file'])


def is_point_inside_single_polygon(pt, polygon, on = True):
    """
    Function that tests if a point is inside a polygon. 
    
    Input:
        pt:         the point to test (a tuple of coordinates)
        polygon:    a ring of the polygon represented by a list of tuple.
        on:         if the point on the boundary is considered inside. 'True' means inside and vice versa.
    Output:
        True:       pt is inside polygon
        False:      pt is outside polygon

    """ 
    flag = False
    for i in range(len(polygon) - 1):
        x = pt[0]
        y = pt[1]
        x1 = polygon[i][0]
        y1 = polygon[i][1]
        x2 = polygon[i + 1][0]
        y2 = polygon[i + 1][1]

        if on and ((x == x1 and y ==y1) or (x == x2 and y ==y2)):  #pt is on the border (one of the vertices) of the polygon
            return True
        
        if y1 <= y < y2 or y2 <= y < y1:
            x_crossing_point = (x2 - x1) / (y2 - y1) * (y - y1) + x1
            if on and x_crossing_point == x:  #pt is on the border (not horizontal lines) of the polygon
                return True                  
            elif x_crossing_point > x:        #pt is inside (not on the border of) the polygon
                flag = not flag              
        elif on and y1 == y == y2:           
            if x1 < x < x2 or x2 < x < x1:    #pt is on the border (horizontal lines) of the polygon
                return True
    return flag


# function describes the experimental variogram
def g(h):
    return 1500 * (1 - math.exp(-9*h**2 / 320**2)) + 0;

def dis(point1, point2):
	return math.sqrt((point2[0]-point1[0])*(point2[0]-point1[0])+(point2[1]-point1[1])*(point2[1]-point1[1]))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
