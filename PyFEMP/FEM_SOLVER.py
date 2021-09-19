# Class version of 1D Solver
from sys import dont_write_bytecode
import numpy as np
import matplotlib as mpl
import warnings


class FEM_Simulation:
    '''
    Object that represents a 1D FEM Simulation.
    '''

    def __init__(self, Element, verbose=False):

        self.Element = Element

        # get initialize data from the element
        self.NoElementDim, \
        self.NoElementNodes, \
        self.ElementDofNames, \
        self.NoElementHistory, \
        self.ElementMaterialNames, \
        self.ElementPostNames = Element.Elmt_Init()

        self.NoElementMaterial = len(self.ElementMaterialNames)
        self.NoNodeDofs = len(self.ElementDofNames)

        # general program variables
        self.verbose = verbose
        self.verbose_system = True
        self.state = 0

        # general discretization variables
        self.time = 0.0                     # current time
        self.dt = 1.0                       # time increment gone frome last time
        self.step = 0                       # current step
        self.lambda_load = 0                # global load multiplier
        self.NoElements = 0                 # number of elements
        self.NoNodes = 0                    # number of nodes
        self.NoDofs = 0                     # number of degrees of freedom
        self.XI = 0                         # nodal coordinates
        self.ELEM = 0                       # element connectivity
        self.h_n = 0                        # previous history field
        self.h_t = 0                        # current history field

        # initialize fields for boundary conditions
        self.NBC = []                       # python list to collect natural boundary conditions before analysis
        self.NBC_Indexes = 0                # vector of indexes to the external load vector where a nbc is present
        self.NBC_Values = 0                 # vector of values to be placed in the external load vector for each nbc index
        self.EBC = []                       # python list to collect essential boundary conditions before analysis
        self.EBC_Indexes = 0                # vector of indexes of constrained degrees of freedom
        self.EBC_Values = 0                 # vector of values for each constrained degree of freedom
        self.NoEquations = 0                # number of all unconstrained dofs

        # element discretization parameter
        self.ElementMaterial = []           # list of material parameter
        self.h_n = 0                        # vector of element history field of t=t   (previous)
        self.h_t = 0                        # vector of element history field of t=t+1 (current)
        self.DI = 0                         # vector of degrees of freedom
        self.R_ext = 0                      # vector of external forces
        

        # make some noise
        print("FEM Solver Instance Created")
        if (self.verbose): print("Simulation dimensions:                    ", self.NoElementDim)
        if (self.verbose): print("Number of element nodes:                  ", self.NoElementNodes)
        if (self.verbose): print("Names of nodal degrees of freedom:        ", self.ElementDofNames)
        if (self.verbose): print("Names of element parameters:              ", self.ElementMaterialNames)
        if (self.verbose): print("Names of available postprocessing fields: ", self.ElementPostNames)




    def Add_Mesh(self, NodesList, ElementConnectivity, verbose=False):
        ''' 
        Add_Mesh(self, NodesList, ElementConnectivity, verbose=False) -> void
        Sets a mesh based on a list of nodes and matrix of element connectivity.

        Input :
        NodeList            -> List of nodal coordinates [... , [x,y], ...]
        ElementConnectivity -> Matrix of nodal indexes per element [... , [n1, n2, n3], ...]
        '''
        # check input
        if (NodesList.ndim == 1):
            no_mesh_no, mesh_dim = len(NodesList), 1
        else:
            no_mesh_no, mesh_dim = NodesList.shape

        no_mesh_el, mesh_no_el = ElementConnectivity.shape

        if (self.verbose or verbose): print('Mesh NoNodes          : ',no_mesh_no)
        if (self.verbose or verbose): print('Mesh Dimension        : ',mesh_dim)
        if (self.verbose or verbose): print('Mesh NoElements       : ',no_mesh_el)
        if (self.verbose or verbose): print('Mesh Nodes per Element: ',mesh_no_el)

        if (self.NoElementDim != mesh_dim): raise NameError('Mesh dimension is not the same as elements.')
        if (self.NoElementNodes != mesh_no_el): raise NameError('Mesh is not compatible to element topology.')

        # process infos
        self.NoElements = no_mesh_el
        self.NoNodes = no_mesh_no
        self.NoDofs = no_mesh_no * self.NoNodeDofs
        if (verbose): print('Mesh Total Dofs       : ',self.NoDofs)

        self.XI   = np.array(NodesList, dtype=np.float64)
        self.ELEM = np.array(ElementConnectivity, dtype=np.uint)

        if (self.verbose or self.verbose): print(' Finite Elemenmt Mesh Read!')
        self.state = 1




    def Add_Material(self, MaterialList, Option=None):
        '''Adds Material parameters as a list [....] for the next element without already specified material. With Option=All, all elements are set with the given list of parameters.'''
        if len(MaterialList) != self.NoElementMaterial:
            print(
                'Error: Number of material parameter does not fit element requirements!')
            print('       Requred parameters are :')
            print(*self.ElementMaterialNames, sep=", ")
            raise NameError('Error processing material parameters!')

        if Option == "All":
            self.ElementMaterial = []
            for i in range(self.NoElements):
                self.ElementMaterial.append(MaterialList)
            if (self.verbose): print(' Material set for All Elements')

        if Option == None:
            self.ElementMaterial.append(MaterialList)
            if (self.verbose): print(' Material set for Element %i' % len(self.ElementMaterial))

    
    
    def Add_EBC(self, NodeSelector, DofSelector, Value):
        '''Sets an essential boundary condition by NodeSelector, DofSelector, Value'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            self.EBC.append([node, AffectedDof, Value])

    
    
    def Add_NBC(self, NodeSelector, DofSelector, Value):
        '''Sets an essential boundary condition by NodeSelector, DofSelector, Value'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            self.NBC.append([node, AffectedDof, Value])



    def SelectDof(self, Input):
        '''Returns a single integer for the dof'''
        if isinstance(Input, int):
            return Input
        if isinstance(Input, str):
            for i, dofname in enumerate(self.ElementDofNames):
                if dofname == Input:
                    return i
        raise NameError('Error ! DOF name not supported by element')
        return 100



    def SelectNodes(self, Input):
        '''
        Returns a list containing the node number that fit the input.
        Input can be: 
            A single node index SelectNodes(0)
            A list of indexes node index SelectNodes([0,1,2])
            A conditional based on the dimension 1D: SelectNodes("x==0") 2D: SelectNodes("x==0 && y==0")
        '''
        Outlist = []

        # if input is a singe integer, check if there is a node for this integer and return it in a list
        if isinstance(Input, int):
            if Input <= self.NoElements:
                Outlist.append(Input)
            else:
                print('Error: %i is not a valid node number!' % Input)
                return

        # if input is a list, we check each entry for being an integer and proceed as before. if the entry is
        # an interger we append it to the output list
        if isinstance(Input, list):
            for i in range(len(Input)):
                if isinstance(Input[i], int):
                    if Input[i] <= self.NoElements:
                        Outlist.append(Input[i])
                    else:
                        print('Error: %i is not a valid node number!' %
                              Input[i])
                        return
                else:
                    print(
                        'Error: ', Input[i], " is not a valid node number! Integer required.")

        # if input is a string, it is supposed to be a conditional
        if isinstance(Input, str):

            # 1D - condition is x only
            if (self.NoElementDim==1):
                conditional = eval("lambda x: "+Input)
                Outlist = np.arange(self.NoNodes)[[conditional(x) for x in self.XI]]
            
            # 2D - condition is x and y
            elif (self.NoElementDim==2):
                conditional = eval("lambda x, y: "+Input)
                Outlist = np.arange(self.NoNodes)[[conditional(x,y) for x, y in self.XI]]

        return Outlist



    def Analysis(self):
        '''Enters into the Analysis phase. At least there must be finite elements and Materials'''
        if self.state < 1:
            self.state_report()
            return
        elif self.NoElements < 1:
            raise NameError('Error! No Elements! Use AddMesh.')
        elif len(self.ElementMaterial) != self.NoElements:
            raise NameError('Error! Not sufficent Material provided! Use AddMaterial.')
        
        # initialize history
        self.h_n = np.zeros(self.NoElements * self.NoElementHistory)
        self.h_t = np.copy(self.h_n)

        # initialize degrees of freedom
        self.DI = np.zeros(self.NoNodes * self.NoNodeDofs)

        # initialize external right hand side
        self.R_ext = np.zeros(self.NoNodes * self.NoNodeDofs)

        # consolidate boundary conditions
        self.EBC_Indexes = np.array([ node*self.NoNodeDofs+dof for node, dof, value in self.EBC], dtype=np.uint)
        self.EBC_Values  = np.array([ value for node, dof, value in self.EBC], dtype=np.float64)
        self.NBC_Indexes = np.array([ node*self.NoNodeDofs+dof for node, dof, value in self.NBC], dtype=np.uint)
        self.NBC_Values  = np.array([ value for node, dof, value in self.NBC], dtype=np.float64)
        self.NoEquations = self.NoNodes * self.NoNodeDofs - len(self.NBC_Indexes)

        print('Entering Analysis phase')
        if (self.verbose_system):
            print('---------------------------------')
            print('FE Setup Summary :')
            print('NoElementNodes   :', self.NoElementNodes)
            print('NoNodeDofs       :', self.NoNodeDofs)
            print('ElementDofNames  :', self.ElementDofNames)
            print('ElementPostNames :', self.ElementPostNames)
            print('NoElementHistory :', self.NoElementHistory)
            print('NoElements       :', self.NoElements)
            print('NoNodes          :', self.NoNodes)
            print('NoDofs           :', self.NoDofs)
            print('NoEssential BC   :', len(self.EBC))
            print('NoNatural BC     :', len(self.NBC))
            print('---------------------------------')

        self.state = 100



    def state_report(self):
        '''Gives hints to the user what to do next, based on a standard procedure.'''
        if self.state == 0:
            print('Input required: Call the AddMesh() function.')
        elif self.state == 1:
            print('state is 1')



    def NextStep(self, time=1, lambda_load=1):

        # check requirements
        if self.state < 100:
            print('Error: Simulation has not entered analysis phase via Analysis().')
        if (self.verbose_system):
            print('\nCurrent Time : %f5.2' % time)
        self.dt = time - self.time
        if self.dt < 1e-16:
            print(
                'Error: Time given in NextStep is smaller than internal time: %f5.2' % time)
            return
        
        # time shift time dependent variables
        self.time = time
        self.step += 1
        self.h_n = np.copy(self.h_t)
        self.lambda_load = lambda_load

        # apply EBC to DI
        self.DI[self.EBC_Indexes] = self.lambda_load * self.EBC_Values

        # apply NBC to 
        self.R_ext[self.NBC_Indexes] = self.lambda_load * self.NBC_Values

        return




    def CallElement(self, i, verbose=False):
        if i > self.NoElements:
            print('Error: Input exceeds number of elements. max is : %i8' %
                  self.NoElements)
        elmt_nodes        = self.ELEM[i]
        elmt_dof_indexes  = np.array([i * self.NoNodeDofs + d for i in elmt_nodes for d in range(self.NoNodeDofs)], dtype=np.uint)
        elmt_hist_indexes = np.arange(i * self.NoElementHistory,(i+1) * self.NoElementHistory)

        Elmt_XI = (self.XI[elmt_nodes]).flatten()
        Elmt_UI = self.DI[elmt_dof_indexes]
        Elmt_Hn = self.h_n[elmt_hist_indexes]
        Elmt_Ht = self.h_t[elmt_hist_indexes]
        Elmt_Mat = self.ElementMaterial[i]
        if (self.verbose or verbose):
            print('Calling  : ', i)
            print('Elmt XI  : ', Elmt_XI)
            print('Elmt UI  : ', Elmt_UI)
            print('Elmt Mat : ', Elmt_Mat)
            print('Elmt Hn  : ', Elmt_Hn)
            print('Elmt Ht  : ', Elmt_Ht)
        # call element routine to get element vector / element matrix
        r_e, K_e = self.Element.Elmt_KS(
            Elmt_XI,
            Elmt_UI,
            Elmt_Hn,
            Elmt_Ht,
            Elmt_Mat,
            self.dt
        )
        self.h_t[elmt_hist_indexes] = Elmt_Ht
        if (self.verbose or verbose):
            print('Element Vector :', r_e)
            print('Element Matrix :')
            print(K_e)
        return r_e, K_e



    def Assemble(self, verbose=False):
        r = np.zeros(self.NoDofs)
        K = np.zeros((self.NoDofs, self.NoDofs))
        for e in range(self.NoElements):
            r_e, K_e = self.CallElement(e)

           # Assemble global vector and global matrix
            # compute dof indexes
            dof_indexes = np.array([ n*self.NoNodeDofs + d for n in self.ELEM[e] for d in range(self.NoNodeDofs)] ,dtype=np.uint)

            # assemble right hand side
            r[dof_indexes] = r[dof_indexes] + r_e

            # assemble stiffnes matrix
            for n, i in enumerate(dof_indexes):
                for m, j in enumerate(dof_indexes):
                    K[i][j] += K_e[n][m]

        if (self.verbose or verbose):
            print('Global Vector :', r)
            print('Global Matrix :')
            print(K)
        return r, K



    def FormLinearSystem(self):
        r, K = self.Assemble()
        # Apply Essential Boundary Conditions / Build reduction operator
        ebc_dof_indexes = self.EBC_Indexes

        #    Build reduction operator
        I_full = np.diag(np.ones(self.NoDofs))
        I_red = np.delete(I_full, ebc_dof_indexes, 0)

        # Build reduced System of equations
        #    Sum up force vectors
        RHS = self.R_ext - r

        # Reduce the System
        RHS = I_red.dot(RHS)
        LHS = I_red.dot(K.dot(I_red.T))
        if (self.verbose):
            print('Lin EqS RHS :', RHS)
            print('Lin EqS LHS :')
            print(LHS)

        return RHS, LHS, I_red 



    def NewtonIteration(self):
        '''
        Performes a NewtonIteration on the current state of the system.
        Returns a convergence indicator: sqrt(|R|*|dDI|)
        '''
        RHS, LHS, I_red = self.FormLinearSystem()
        Residual = np.sqrt(RHS.dot(RHS))
        if (self.verbose_system):
            print('      |R|    : %f10.8' % Residual)

        # Solve the linear System
        dDI = np.linalg.solve(LHS, RHS)
        if (self.verbose):
            print('Lin EqS Sol :', dDI)

        # Compute Norm of Solution Vector
        Norm_dDI = np.sqrt(dDI.dot(dDI))
        if (self.verbose_system):
            print('      |dDI|  : %f10.8' % Norm_dDI)

        # Update Global Vector of DOFs
        self.DI += dDI.dot(I_red)
        if (self.verbose):
            print('Current vector of unknowns :', self.DI)

        return np.sqrt(Residual*Norm_dDI)



    def Mod_Material(self, MaterialList, Element_I):
        ''' Modifies an already defined Material of element I'''
        if len(MaterialList) != self.NoElementMaterial:
            print(
                'Error: Number of material parameter does not fit element requirements!')
            print('       Requred parameters are :')
            print(*self.ElementMaterialNames, sep=", ")
            return
        if Element_I > len(self.ElementMaterial) - 1:
            print('Error: Material of Element %i is not defined yet.' % Element_I)
            return

        self.ElementMaterial[Element_I] = MaterialList
        return



    def NodalDof(self, NodeSelector, DofSelector):
        '''Return current DoF by NodeSelector, DofSelector'''
        NodeList = self.SelectNodes(NodeSelector)
        AffectedDof = self.SelectDof(DofSelector)
        if AffectedDof >= self.NoNodeDofs:
            print("Error: Nodal degrees of freedom do not exceed %i" %
                  self.NoNodeDofs)
        for node in NodeList:
            dof_index = (node * self.NoNodeDofs) + AffectedDof
            return self.DI[dof_index]



    def CallElementPost(self, i, PostName):
        '''
        Calls the postprocessing routine Elmt_Post for element i with the current Simulation fields.
        Returns first a list of all elment node indexes and next the vector with the requested PostName
        data, one scalar for each node.
        '''
        if i > self.NoElements:
            print('Error: Input exceeds number of elements. max is : %i8' %
                  self.NoElements)
        elmt_nodes        = self.ELEM[i]
        elmt_dof_indexes  = np.array([i * self.NoNodeDofs + d for i in elmt_nodes for d in range(self.NoNodeDofs)], dtype=np.uint)
        elmt_hist_indexes = np.arange(i * self.NoElementHistory,(i+1) * self.NoElementHistory)

        Elmt_XI = (self.XI[elmt_nodes]).flatten()
        Elmt_UI = self.DI[elmt_dof_indexes]
        Elmt_Hn = self.h_n[elmt_hist_indexes]
        Elmt_Ht = self.h_t[elmt_hist_indexes]
        Elmt_Mat = self.ElementMaterial[i]

        # call element routine to get element vector / element matrix
        r_post_e= self.Element.Elmt_Post(
            Elmt_XI,
            Elmt_UI,
            Elmt_Hn,
            Elmt_Ht,
            Elmt_Mat,
            self.dt,
            PostName
        )
        return elmt_nodes, r_post_e



    def Assemble_Post(self, PostName):
        post_vector = np.zeros(self.NoNodes)

        for e in range(self.NoElements):
            elmt_nodes, r_post_e = self.CallElementPost(e, PostName)
            
            # post weights
            post_weights = 1/np.array([np.count_nonzero(self.ELEM== i) for i in elmt_nodes])

            # assemble post vector
            post_vector[elmt_nodes] = post_vector[elmt_nodes] + (r_post_e * post_weights)

        return post_vector



    def PostProcessing(self, PostName, points = []):
        '''
        Returns a list of nodes and a list of 
        one requested scalar for each of these nodes.
        The requested scalar is computed from the element 
        subroutine "Elmt_Post" by PostName.
        By default, the PostName field is returned for all mesh nodes.
        x -> matrix that contains the nodal positions
        p -> vector of values for each node
        Example:
            PostProcessing("UX") -> x, p
        Optionally, specify arbitrary positions for which to evaluate the PostNames
        Examples:
            PostProcessing("UX", [0.0, 0.0]) -> x, p
            PostProcessing("UX", [[0.0, 0.0],...,[1.0, 0.0]]) -> x, p
        '''
        import matplotlib.tri as mtri
        if PostName not in self.ElementPostNames: 
            print('Warning, PostName not available.')
            print('Choose from: ') 
            print( [str(name)+", " for name in self.ElementPostNames])
            return
        if (len(points)==0):
            post_vector = self.Assemble_Post(PostName)
            return self.XI.copy(), post_vector
        else:
            # assume points is an input of points
            if (self.NoElementDim==2):
                post_vector = self.Assemble_Post(PostName)
                if (self.NoElementNodes==3):
                    XI = self.XI
                    Elmt = self.ELEM
                elif (self.NoElementNodes==4):
                    XI = self.XI
                    Elmt = np.array([[elmt[0], elmt[1], elmt[3], elmt[1], elmt[2], elmt[3]] for elmt in self.ELEM] ,dtype=np.uint)
                    Elmt = Elmt.reshape(-1,3)
                else:
                    raise NameError('Error! No 2D postprocessing for :'+str(self.NoElementNodes)+' nodes')

                mesh = mtri.Triangulation(XI[:,0], XI[:,1], Elmt)
                interpolated_data = mtri.LinearTriInterpolator(mesh, post_vector)

                if (len(np.shape(points))==1): # single point
                    post_vector = interpolated_data(points[0], points[1])
                else:
                    post_vector = interpolated_data(points[:,0], points[:,1])

                return points, post_vector.filled(0.0)

                
            else:
                raise NameError('Error! No xD postprocessing for :'+str(self.NoElementNodes)+' nodes')
                return 0, 0



    def ShowMesh(self, ax, deformedmesh = False, boundaryconditions = False, PostName = "", **kwargs):
        '''
        Visualisation of Mesh and solution.
        ShowMesh(self, ax, deformedmesh = False, boundaryconditions = False, PostName = "", **kwargs) -> plot
        The function takes a matplotlit figure axis as arguments and plots inside the same.
        The option deformedmesh = False can be used to plot the deformed mesh, where
        nodal displacement is supposed to be the n.th first nodal degrees of freedom
        in n-dimensions.
        PostName = "" can be used to create a countour plot of the quantities provided by the
        element subroutine "Elmt_Post".  **kwargs is passed to the sub-plotting routines
        (patches.Polygon and tricontourf) to modify eg line colors/witdhs etc.
        '''
        import matplotlib as mpl

        if (self.NoElementDim==2):
            # 2D Visualisation

            postplot = 0
            
            if (self.NoElementNodes==3):
                #  visualisation of T1
                if self.verbose: print('Visualisation 2D T1')
                XI = np.copy(self.XI)
                
                # some customization on optional arguments
                if 'ec' not in kwargs.keys(): kwargs['ec']=(0.2, 0.2, 0.2, 1.0)
                if 'label' in kwargs.keys():
                    ax.scatter([],[],label=kwargs['label'], ec=kwargs['ec'], fc=kwargs['ec'])
                    del kwargs['label']

                
                if deformedmesh:
                    Ux = self.DI[::self.NoNodeDofs]
                    Uy = self.DI[1::self.NoNodeDofs]
                    XI = XI + np.array([Ux, Uy]).T
                
                for elem in self.ELEM:
                    ax.add_patch(mpl.patches.Polygon(
                        XI[elem],
                        True,
                        fc=(0.0, 0.0, 0.0, 0.0),
                        **kwargs
                        ))
                
                ax.set_xlim(min(XI[:,0])*1.1, max(XI[:,0])*1.1)
                ax.set_ylim(min(XI[:,1]*1.1), max(XI[:,1])*1.1)
                ax.set_aspect('equal')

                if (PostName!=""):
                    if PostName not in self.ElementPostNames: 
                        print('Warning, PostName not available.')
                        print('Choose from: ') 
                        print( [str(name)+", " for name in self.ElementPostNames])
                        return

                    # get the post vector
                    post_vector = self.Assemble_Post(PostName)
                    
                    # counter plot
                    warnings.filterwarnings("ignore") # to supress a warning from countour
                    levels = 5
                    postplot = ax.tricontourf(XI[:,0], XI[:,1], post_vector, levels, triangles=self.ELEM, **kwargs, cmap='jet')
                    warnings.filterwarnings("default")

                

            elif (self.NoElementNodes==4):
                
                #  visualisation of Q1
                if self.verbose: print('Visualisation 2D Q1')
                XI = np.copy(self.XI)
                
                # some customization on optional arguments
                if 'ec' not in kwargs.keys(): kwargs['ec']=(0.2, 0.2, 0.2, 1.0)
                if 'label' in kwargs.keys():
                    ax.scatter([],[],label=kwargs['label'], ec=kwargs['ec'], fc=kwargs['ec'])
                    del kwargs['label']


                if deformedmesh:
                    Ux = self.DI[::self.NoNodeDofs]
                    Uy = self.DI[1::self.NoNodeDofs]
                    XI = XI + np.array([Ux, Uy]).T

                for elem in self.ELEM:
                    ax.add_patch(mpl.patches.Polygon(
                        XI[elem],
                        True,
                        fc=(0.0, 0.0, 0.0, 0.0),
                        **kwargs
                        ))

                ax.set_xlim(min(XI[:,0])*1.1, max(XI[:,0])*1.1)
                ax.set_ylim(min(XI[:,1]*1.1), max(XI[:,1])*1.1)
                ax.set_aspect('equal')

                if (PostName!=""):
                    if PostName not in self.ElementPostNames: 
                        print('Warning, PostName not available.')
                        print('Choose from: ') 
                        print( [str(name)+", " for name in self.ElementPostNames])
                        return

                    # get the post vector
                    post_vector = self.Assemble_Post(PostName)
                    
                    # for Q1 topology we build a triangulation
                    triangulation_for_Q1 = np.array([[elmt[0], elmt[1], elmt[3], elmt[1], elmt[2], elmt[3]] for elmt in self.ELEM] ,dtype=np.uint)
                    triangulation_for_Q1 = triangulation_for_Q1.reshape(-1,3)
                    
                    # counter plot
                    warnings.filterwarnings("ignore") # to supress a warning from countour
                    levels = 5
                    postplot = ax.tricontourf(XI[:,0], XI[:,1], post_vector, levels, triangles=triangulation_for_Q1, **kwargs)
                    warnings.filterwarnings("default")

            else:
                raise NameError('Error! No 2D visualisation for :'+str(self.NoElementNodes)+' nodes')

            if (boundaryconditions):
                dofs_per_node = self.NoNodeDofs
                if (dofs_per_node>7): 
                    raise NameError("Error! Visualisation of boundaray conditions for so many dofs not supported!")

                class bc_vis:
                    def __init__(self):
                        self.x = []
                        self.y = []

                ebc_label=['EBC DOF 0', 'EBC DOF 1', 'EBC DOF 2', 'EBC DOF 3', 'EBC DOF 4', 'EBC DOF 5', 'EBC DOF 6']
                ebc_marker=['v', '^', 'D', 's', 'p', '+', '*']
                ebc_visualisations = [bc_vis() for i in range(dofs_per_node)]
                for ebc in self.EBC:
                    affected_node = ebc[0]
                    affected_dof = ebc[1]
                    affected_XI = XI[affected_node] # is provided by steps before also in deformed version
                    ebc_visualisations[affected_dof].x.append(affected_XI[0])
                    ebc_visualisations[affected_dof].y.append(affected_XI[1])

                nbc_label=['NBC DOF 0', 'NBC DOF 1', 'NBC DOF 2', 'NBC DOF 3', 'NBC DOF 4', 'NBC DOF 5', 'NBC DOF 6']
                nbc_marker=["4", "2", 'D', 's', 'p', '+', '*']
                nbc_visualisations = [bc_vis() for i in range(dofs_per_node)]
                for nbc in self.NBC:
                    affected_node = nbc[0]
                    affected_dof = nbc[1]
                    affected_XI = XI[affected_node]
                    nbc_visualisations[affected_dof].x.append(affected_XI[0])
                    nbc_visualisations[affected_dof].y.append(affected_XI[1])

                for i in range(dofs_per_node):
                    ax.plot(
                        ebc_visualisations[i].x, 
                        ebc_visualisations[i].y, 
                        markeredgecolor='red', marker=ebc_marker[i], 
                        linestyle='None', markerfacecolor='None', label=ebc_label[i])
                for i in range(dofs_per_node):
                    ax.plot(
                        nbc_visualisations[i].x, 
                        nbc_visualisations[i].y, 
                        markeredgecolor='blue', marker=nbc_marker[i], 
                        linestyle='None', markerfacecolor='None', label=nbc_label[i])

            return postplot


        else: 
            raise NameError("Error! Visualisation not supported!")
        return
