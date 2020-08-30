def msh_rec(X0, X1, N, type='Q1'):
    """
    msh_rec(X0, X1, N, type='Q1') -> XI, ELEM
    Function to return a regular mesh from a rectangle.
    The mesh either consists of linear rectangle elements with 4 nodes (type='Q1' (default))
    or linear triangles with 3nodes (type='T1').

    Input:
    X0           -> 2D coordinates of lower left corner
    X1           -> 2D coordinates of upper right corner
    N            -> regular mesh divisioning [nx, ny]
    type = 'Q1' -> specify mesh topology

    Output:
    XI           -> np.array with the nodal coordinates [...,[x,y],...]
    ELEM         -> np.array with the nodal indexes for each element [...,[n1,n2,n3,n4],...]
    """
    import numpy as np

    nx, ny = N

    # compute length
    dx = (X1[0] - X0[0]) / nx
    dy = (X1[1] - X0[1]) / ny

    # compute nodal coordinates
    x = np.linspace(X0[0], X1[0], nx+1, dtype=np.float64)
    y = np.linspace(X0[1], X1[1], ny+1, dtype=np.float64)
    XI = np.array([np.array([xx, yy], dtype=np.float64) for yy in y for xx in x], dtype=np.float64)

    # generate elements
    if (type=='Q1'):
        ELEM = np.array([np.array([ \
            nyi     *(nx+1)+nxi  ,  \
            nyi     *(nx+1)+nxi+1,  \
            (nyi+1) *(nx+1)+nxi+1,  \
            (nyi+1) *(nx+1)+nxi     \
                 ], dtype=np.uint) for nyi in range(ny) for nxi in range(nx)], dtype=np.uint)

    elif (type=='T1'):
        ELEM = np.array([np.array([ \
            nyi     *(nx+1)+nxi  ,  \
            nyi     *(nx+1)+nxi+1,  \
            (nyi+1) *(nx+1)+nxi,    \
            nyi     *(nx+1)+nxi+1,  \
            (nyi+1) *(nx+1)+nxi+1,  \
            (nyi+1) *(nx+1)+nxi     \
                 ], dtype=np.uint) for nyi in range(ny) for nxi in range(nx)], dtype=np.uint).reshape((-1,3))
    else:
        raise NameError('Type :'+type+" not supported!")
        return 0, 0

    return XI, ELEM



def msh_conv_quad(X1, X2, X3, X4, N, type='Q1'):
    """
    msh_conv_quad(X1, X2, X3, X4, N, type='Q1') -> XI, ELEM
    Function to return a regular mesh from a 4 noded polygon.
    The mesh either consists of linear rectangle elements with 4 nodes (type='Q1' (default))
    or linear triangles with 3nodes (type='T1').

    Input:
    X0           -> 2D coordinates of lower left corner
    X1           -> 2D coordinates of lower right corner
    X2           -> 2D coordinates of upper right corner
    X3           -> 2D coordinates of upper left corner
    N            -> regular mesh divisioning [nx, ny]
    type = 'Q1' -> specify mesh topology

    Output:
    XI           -> np.array with the nodal coordinates [...,[x,y],...]
    ELEM         -> np.array with the nodal indexes for each element [...,[n1,n2,n3,n4],...]
    """
    import numpy as np

    
    # get undeformed mesh
    UXI, ELEM = msh_rec([0.0, 0.0], [1.0, 1.0], [N[0], N[1]], type=type)

    # map each nodal coordinate on the actual poygon by linear interpolation
    XI = np.array([np.array([
        ((1-x))*((1-y))*X1[0] +(x)*((1-y))*X2[0] + ((1-x))*(y)*X4[0] + (x)*(y)*X3[0],\
        ((1-x))*((1-y))*X1[1] +(x)*((1-y))*X2[1] + ((1-x))*(y)*X4[1] + (x)*(y)*X3[1]\
    ], dtype=np.float64) for x, y in UXI], dtype=np.float64)


    return XI, ELEM





# Test:
if __name__ == "__main__":

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    XI, ELEM = msh_conv_quad([0.0, 0.0], [48.0, 44.0], [48.0, 60.0], [0.0, 44.0], [10, 4], type='T1')

    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    ax.scatter(XI[:,0], XI[:,1])
    for element in ELEM:
        ax.add_patch(mpl.patches.Polygon(XI[element], True, fc=(0.0, 0.0, 0.0, 0.0), ec='black'))
    plt.show()

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    XI, ELEM = msh_rec([0.0, 0.0], [10.0, 5.0], [10, 5], type='T1')

    fig, ax = plt.subplots(1, 1, figsize = (10, 5))
    ax.scatter(XI[:,0], XI[:,1])
    for element in ELEM:
        ax.add_patch(mpl.patches.Polygon(XI[element], True, fc=(0.0, 0.0, 0.0, 0.0), ec='black'))
    plt.show()