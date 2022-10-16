from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import copy as cp
# np.randomz.seed(12)

class Node:
    def __init__(self,l,r):
        self.l_nei = None
        self.r_nei = None

        self.l = l
        self.r = r
    
    def get_loop(self,edges):
        boundary = np.array([self.r])
        if self.r_nei != None:
            boundary = np.concatenate((boundary, edges[self.r_nei].get_loop(edges)),axis=0)
        else:
            return boundary
        return boundary

class Grid:
    def __init__(self,tl,bl,tr,br):
        self.l_nei = None
        self.t_nei = None
        self.r_nei = None
        self.b_nei = None

        self.tl = tl
        self.bl = bl
        self.tr = tr
        self.br = br
        
        self.centroid = (self.tl+self.bl+self.tr+self.br)/4.0
    
    def get_centroid(self):
        self.centroid = (self.tl+self.bl+self.tr+self.br)/4.0

    def get_x(self):
        return np.array([self.bl[0],self.br[0],self.tr[0],self.tl[0],self.bl[0]])

    def get_y(self):
        return np.array([self.bl[1],self.br[1],self.tr[1],self.tl[1],self.bl[1]])

    def get_corners(self):
        return np.array([self.bl,self.br,self.tr,self.tl])

# class Node:
#     def __init__(self,id,arena_w,arena_h):
#         self.node_id = np.array([int(id/arena_w),id%arena_w])

def adjacency_mat(grids):
    inner_mat = np.zeros((grids.shape[1],grids.shape[1]))
    outter_mat = np.zeros((grids.shape[0]*grids.shape[1],grids.shape[0]*grids.shape[1]))

    for i in range(grids.shape[1]):
        if i+1<grids.shape[1]:
            inner_mat[i,i+1] = 1
        if i-1>=0:
            inner_mat[i,i-1] = 1
    
    for i in range(outter_mat.shape[1]):
        if i+grids.shape[1]<outter_mat.shape[1]:
            outter_mat[i,i+grids.shape[1]] = 1
        if i-grids.shape[1]>=0:
            outter_mat[i,i-grids.shape[1]] = 1
    
    vertices = np.kron(np.eye(grids.shape[0]),inner_mat) + outter_mat
    return vertices

def get_plot(array,type_):
    if type_ == graph:
        x = np.empty(0)
        y = np.empty(0)
        for i in range(len(array)):
            x,y = np.append(x,array[i][0]),np.append(y,array[i][1])
            plt.plot(x,y)
            
    elif type_ == edge:
        count = 0
        x = np.empty(0)
        y = np.empty(0)
        for elem in array:
            x,y = np.append(x,elem[0]),np.append(y,elem[1])
            count+=1
            if count == 2:
                plt.plot(x,y)
                x = np.empty(0)
                y = np.empty(0)
                count = 0

def plot_show():
    plt.figure(figsize = (5,5))
    plt.show()

def get_coordinates(vertices,boundry_check,action,lim_x,lim_y,edge_length):

    cur_x = vertices[-1][0]
    cur_y = vertices[-1][1]

    boundry = np.array([0,0,0,0])

    if(cur_x == lim_x):
        boundry[right] = 1
        boundry_check[right] = 1
    if(cur_x == 0):
        boundry[left] = 1
        boundry_check[left] = 1
    if(cur_y == lim_y):
        boundry[up] = 1
        boundry_check[up] = 1
    if(cur_y == 0):
        boundry[down] = 1
        boundry_check[down] = 1

    arr =  (boundry_check == 1)
    boundry_condition = True

    for value in arr:
        if value == False:
            boundry_condition = False

    if action == up:
        next_x = cur_x 
        next_y = cur_y+1
        collision = False
        for i in range(len(vertices)-1):
            if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
                collision=True

        if(boundry[up] == True):
            return False,edge_length,vertices

        elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            edge_length+=1

            return True,edge_length,vertices

        elif collision:
            return False,edge_length,vertices

        else:
            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            updated_edge_len = edge_length+1

            trial = 0
            flag = True
            actions = [0,1,2,3]
            while flag:
                trial +=1
                actions = [0,1,2,3]
                choose_action = sample(actions,1)
                actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
                fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

                if fact_vect_temp == True:
                    indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
                    unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
                    if np.array_equal(unique_vertices,temp_vertices[1:]):
                        return True,edge_len_temp,temp_vertices
                if trial == 4:
                    return False,edge_length,vertices
                else:
                    flag = True

    elif action == right:
        next_x = cur_x+1 
        next_y = cur_y

        collision = False
        for i in range(len(vertices)-1):
            if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
                collision=True

        if(boundry[right] == True):

            return False,edge_length,vertices

        elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10  ):

            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            edge_length+=1

            return True,edge_length,vertices

        elif(collision):
            return False,edge_length,vertices

        else:
            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            updated_edge_len = edge_length+1

            trial = 0
            flag = True
            actions = [0,1,2,3]
            while flag:
                trial +=1
                actions = [0,1,2,3]
                choose_action = sample(actions,1)
                actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
                fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

                if fact_vect_temp == True:
                    indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
                    unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
                    if np.array_equal(unique_vertices,temp_vertices[1:]):
                        return True,edge_len_temp,temp_vertices
                if trial == 4:
                    return False,edge_length,vertices
                else:
                    flag = True

    elif action == down:
        next_x = cur_x 
        next_y = cur_y-1
        collision = False
        for i in range(len(vertices)-1):
            if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
                collision=True

        if(boundry[down] == True):
            return False,edge_length,vertices

        elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            edge_length+=1

            return True,edge_length,vertices

        elif(collision):
            return False,edge_length,vertices

        else:
            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            updated_edge_len=edge_length+1                  

            trial = 0
            flag = True
            actions = [0,1,2,3]
            while flag:
                trial +=1
                
                choose_action = sample(actions,1)
                actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
                fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

                if fact_vect_temp == True:
                    indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
                    unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
                    if np.array_equal(unique_vertices,temp_vertices[1:]):
                        return True,edge_len_temp,temp_vertices
                if trial == 4:
                    return False,edge_length,vertices
                else:
                    flag = True

    elif action == left:
        next_x = cur_x-1 
        next_y = cur_y
        collision = False
        for i in range(len(vertices)-1):
            if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
                collision=True

        if(boundry[left] == True):
            return False,edge_length,vertices

        elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

            vertices = np.append(vertices,[[next_x,next_y]],axis=0)
            edge_length+=1

            return True,edge_length,vertices

        elif(collision):
            return False,edge_length,vertices

        else:
            vertices_ = np.append(vertices,[[next_x,next_y]],axis=0)
            updated_edge_len=edge_length+1                 

            trial = 0
            flag = True
            actions = [0,1,2,3]
            while flag:
                trial +=1
                
                choose_action = sample(actions,1)
                actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
                fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices_,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

                if fact_vect_temp == True:
                    indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
                    unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
                    if np.array_equal(unique_vertices,temp_vertices[1:]):
                        return True,edge_len_temp,temp_vertices
                if trial == 4:
                    return False,edge_length,vertices
                else:
                    flag = True

def region_generator(length,breadth):
    # np.random.seed()
    nodes = []
    t_edges = []
    b_edges = []
    l_edges = []
    r_edges = []
    for b in range(0,breadth+1,int(breadth/10)):
        for l in range(0,length+1,int(length/10)):
            nodes.append([b,l])
            if b == 0:
                t_edges.append([l,breadth])
                b_edges.append([l,0])        
        r_edges.append([length,b])
        l_edges.append([0,b])

    edges = np.concatenate((np.array(b_edges),np.array(r_edges[1:]),np.flip(np.array(t_edges),axis=0)[1:],np.flip(np.array(l_edges),axis=0)[1:-1]),axis=0)
    fig,ax = plt.subplots()
    ax.set_aspect('equal')
    ax.plot(edges[:,0],edges[:,1],color='red')
    ax.plot([edges[0,0] for i in range(2)],np.array(l_edges)[0:2,1],color='red')
    ax.scatter(edges[:,0],edges[:,1])
    plt.ion()
    plt.show()
    times = np.random.randint(0,100)
    
    for run in range(times):
        delete_edge_loc = np.random.randint(0,len(edges)-1)
        deleted_edge = np.copy(edges[delete_edge_loc])
        ax.plot([edges[i,0] for i in range(delete_edge_loc-1,delete_edge_loc+2)],[edges[i,1] for i in range(delete_edge_loc-1,delete_edge_loc+2)],color='black')
        plt.show()
        plt.pause(0.1)
        append_edges_loc = delete_edge_loc
        ax.scatter(deleted_edge[0],deleted_edge[1],c='black')
        plt.show()
        plt.pause(0.1)
        vec_1 = edges[delete_edge_loc-1] - edges[delete_edge_loc]
        ax.plot([edges[delete_edge_loc-1][0], edges[delete_edge_loc][0]],[edges[delete_edge_loc-1][1], edges[delete_edge_loc][1]],color='orange')
        vec_2 = edges[delete_edge_loc+1] - edges[delete_edge_loc]
        ax.plot([edges[delete_edge_loc+1][0], edges[delete_edge_loc][0]],[edges[delete_edge_loc+1][1], edges[delete_edge_loc][1]],color='orange')
        dot = np.round(np.dot(vec_1,vec_2))
        check_line = edges[delete_edge_loc-1] - edges[delete_edge_loc+1]
        cross = np.cross(vec_2,check_line)
        
        if dot==0 and cross>=0: # Cross>0, dot=0
            edges = np.delete(edges,delete_edge_loc,axis=0)
            if (b_edges == deleted_edge).all(1).any() and delete_edge_loc!=0:
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc]+np.array([-length/10,0])],axis=0)
            elif (r_edges == deleted_edge).all(1).any():
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc]+np.array([0,-breadth/10])],axis=0)
            elif (t_edges == deleted_edge).all(1).any():
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc]+np.array([length/10,0])],axis=0)
            elif delete_edge_loc==0:
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc]+np.array([0,breadth/10])],axis=0)

            
            ax.clear()
            
            ax.plot(edges[:,0],edges[:,1],color='blue')
            ax.plot([edges[-1,0],edges[0,0]],[edges[-1,1],edges[0,1]],color='blue')
            ax.scatter(edges[:,0],edges[:,1],color='blue')
            ax.scatter([deleted_edge[0]],[deleted_edge[1]],c='black')
            plt.show()
        elif dot==0 and cross<0:
            pass
            # ax.scatter([deleted_edge[0]],[deleted_edge[1]],c='black')
            # edges = np.delete(edges,delete_edge_loc-1,axis=0)

            # edges = np.delete(edges,delete_edge_loc-1,axis=0)
            # edges = np.insert(edges,append_edges_loc-1,[edges[delete_edge_loc-2]+np.array([-breadth/10,0]),edges[delete_edge_loc-2]+np.array([-2*breadth/10,0]),edges[delete_edge_loc-2]+np.array([-2*breadth/10,length/10])],axis=0)
            # ax.clear()
            # ax.plot(edges[:,0],edges[:,1])
            # ax.scatter(edges[:,0],edges[:,1])
            
            # plt.show()
        elif dot>0:# cross 0, dot >0 ; cross 0, dot <0 
            edges = np.delete(edges,delete_edge_loc,axis=0)

            ax.clear()
            # plt.ion()
            # for i in range(len(edges[:,0])-1):
            #     ax.plot([edges[i,0],edges[i+1,0]],[edges[i,1],edges[i+1,1]],color='orange')
            #     plt.show()
            #     plt.pause(0.5)
            # plt.ioff()
            ax.plot(edges[:,0],edges[:,1],color='blue')
            ax.plot([edges[-1,0],edges[0,0]],[edges[-1,1],edges[0,1]],color='blue')
            ax.scatter(edges[:,0],edges[:,1],color='blue')
            ax.scatter([deleted_edge[0]],[deleted_edge[1]],c='black')
            plt.show()
        elif dot<0:# cross 0, dot >0 ; cross 0, dot <0 
            edges = np.delete(edges,delete_edge_loc,axis=0)
            add_sub = np.random.choice([-1,1])
            if (t_edges == deleted_edge).all(1).any() or (b_edges == deleted_edge).all(1).any():
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc-1]+np.array([0,add_sub*int(breadth/10)]),deleted_edge+np.array([0,add_sub*int(breadth/10)]),edges[delete_edge_loc]+np.array([0,add_sub*int(breadth/10)])],axis=0)
            elif (l_edges == deleted_edge).all(1).any() or (r_edges == deleted_edge).all(1).any():
                edges = np.insert(edges,append_edges_loc,[edges[delete_edge_loc-1]+np.array([add_sub*int(breadth/10),0]),deleted_edge+np.array([add_sub*int(breadth/10),0]),edges[delete_edge_loc]+np.array([add_sub*int(breadth/10),0])],axis=0)
            
            ax.clear()
            # plt.ion()
            # for i in range(len(edges[:,0])-1):
            #     ax.plot([edges[i,0],edges[i+1,0]],[edges[i,1],edges[i+1,1]],color='orange')
            #     plt.show()
            #     plt.pause(0.5)
            # plt.ioff()
            ax.plot(edges[:,0],edges[:,1],color='blue')
            ax.plot([edges[-1,0],edges[0,0]],[edges[-1,1],edges[0,1]],color='blue')
            ax.scatter(edges[:,0],edges[:,1],color='blue')
            ax.scatter([deleted_edge[0]],[deleted_edge[1]],c='black')
            plt.show()

    deletion_list = []
    for i in range(len(edges)-2,1,-2):
        vec_1 = edges[i-1] - edges[i]
        ax.plot([edges[i-1][0], edges[i][0]],[edges[i-1][1], edges[i][1]],color='orange')
        vec_2 = edges[i+1] - edges[i]
        ax.plot([edges[i+1][0], edges[i][0]],[edges[i+1][1], edges[i][1]],color='orange')
        dot = np.round(np.dot(vec_1,vec_2))
        cross = np.cross(vec_2,vec_1)
        if np.linalg.norm(edges[i-1]-edges[i+1]) == 0:
            deletion = 1
        else:
            deletion = 0
        if dot>0 and cross==0:
            deletion_list.append(i)
        if deletion:
            deletion_list.append(i+1)
    deletion_list = np.unique(deletion_list)
    for j in range(len(deletion_list)-1,0,-1):
        edges = np.delete(edges,deletion_list[j],axis=0)
        ax.clear()
        ax.plot(edges[:,0],edges[:,1],color='blue')
        ax.plot([edges[-1,0],edges[0,0]],[edges[-1,1],edges[0,1]],color='blue')
        ax.scatter(edges[:,0],edges[:,1],color='blue')
        ax.scatter([deleted_edge[0]],[deleted_edge[1]],c='black')
        plt.show()
        plt.pause(0.5)    

# region_generator(500,600)

def Inflate_Cut_algorithm(num_vertices):
    '''
    Input: Num_vertices (Must be even and >=4)
    Returns: Array of vertices forming the random simple rectilinear region
    '''

    def polygon_boundary(cells,ax):
        def get_edge(cell,side,one_vertex = 1):
            if side == 't':
                if one_vertex:
                    return [list(cell.tl)]
                else:
                    return [list(cell.tr),list(cell.tl)]
            elif side == 'b':
                if one_vertex:
                    return [list(cell.br)]
                else:
                    return [list(cell.bl),list(cell.br)]
            elif side == 'r':
                if one_vertex:
                    return [list(cell.tr)]
                else:
                    return [list(cell.br),list(cell.tr)]
            elif side == 'l':
                if one_vertex:
                    return [list(cell.bl)]
                else:
                    return [list(cell.tl),list(cell.bl)]

        grids = cp.copy(cells)
        if len(grids) == 1:
            boundary = None
            for g in grids:
                boundary = grids[g].get_corners()
            boundary = np.insert(boundary,len(boundary),boundary[0],axis=0)
            return boundary
        delete_these = []

        for i in grids:
            marker = [0,0,0,0]  # l,r,t,b
            if grids[i].l_nei is not None:
                marker[0] = 1
            if grids[i].r_nei is not None:
                marker[1] = 1
            if grids[i].t_nei is not None:
                marker[2] = 1
            if grids[i].b_nei is not None:
                marker[3] = 1
            if np.sum(marker) == 4:
                if grids[grids[i].t_nei].l_nei is not None and grids[grids[i].t_nei].r_nei is not None and grids[grids[i].b_nei].l_nei is not None and grids[grids[i].b_nei].r_nei is not None:
                    delete_these.append(i)
        for i in delete_these:
            del grids[i]
        
        # for g in grids:
        #     x = grids[g].get_x()
        #     y = grids[g].get_y()
        #     ax.plot(x,y,color='maroon')
        #     plt.show()

        boundary = []
        edges = {}
        g = list(grids.items())[0][0]
        done = []
        start_edge = None
        end_edge = None
        backup = []
        while True:
            # x = grids[g].get_x()
            # y = grids[g].get_y()
            # ax.plot(x,y,color='yellow')
            # ax.scatter(grids[g].centroid[0],grids[g].centroid[1])
            # plt.show()

            marker = [0,0,0,0]  # ['r','t','l','b']
            if grids[g].l_nei is not None:
                marker[2] = 1
            if grids[g].r_nei is not None:
                marker[0] = 1
            if grids[g].t_nei is not None:
                marker[1] = 1
            if grids[g].b_nei is not None:
                marker[3] = 1
            sides = ['r','t','l','b']

            if len(boundary)==0 or (g != list(grids.items())[0][0] and g not in done):
                
                for s in range(4):
                    if marker[s]==0:
                        edge = get_edge(grids[g],sides[s],one_vertex=0)
                        if len(edges) == 0:
                            start_edge = str(len(edges))
                            end_edge = str(len(edges))
                            edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                        else:
                            for e in [start_edge,end_edge]:
                                ll = np.linalg.norm(np.array(edge[0])-edges[e].l)
                                lr = np.linalg.norm(np.array(edge[0])-edges[e].r)
                                rl = np.linalg.norm(np.array(edge[1])-edges[e].l)
                                rr = np.linalg.norm(np.array(edge[1])-edges[e].r)
                                if ll == 0:
                                    start_edge = str(len(edges))
                                    edges[e].l_nei = str(len(edges))
                                    edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                                    edges[str(len(edges)-1)].r_nei = e
                                    break
                                elif lr == 0:
                                    end_edge = str(len(edges))
                                    edges[e].r_nei = str(len(edges))
                                    edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                                    edges[str(len(edges)-1)].l_nei = e
                                    break
                                elif rl == 0:
                                    start_edge = str(len(edges))
                                    edges[e].l_nei = str(len(edges))
                                    edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                                    edges[str(len(edges)-1)].r_nei = e
                                    break
                                elif rr == 0:
                                    end_edge = str(len(edges))
                                    edges[e].r_nei = str(len(edges))
                                    edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                                    edges[str(len(edges)-1)].l_nei = e
                                    break
                                elif e == end_edge and ll!=0 and rl!=0 and lr!=0 and rr!=0:
                                    backup.append(edge)
                                    break

                if len(backup)!=0:
                    for i in range(len(backup)-1,-1,-1):
                        edge = backup[i]
                        for e in [start_edge,end_edge]:
                            ll = np.linalg.norm(np.array(edge[0])-edges[e].l)
                            lr = np.linalg.norm(np.array(edge[0])-edges[e].r)
                            rl = np.linalg.norm(np.array(edge[1])-edges[e].l)
                            rr = np.linalg.norm(np.array(edge[1])-edges[e].r)
                            if ll == 0:
                                start_edge = str(len(edges))
                                edges[e].l_nei = str(len(edges))
                                edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                                edges[str(len(edges)-1)].r_nei = e
                                del backup[i]
                                break
                            elif lr == 0:
                                end_edge = str(len(edges))
                                edges[e].r_nei = str(len(edges))
                                edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                                edges[str(len(edges)-1)].l_nei = e
                                del backup[i]
                                break
                            elif rl == 0:
                                start_edge = str(len(edges))
                                edges[e].l_nei = str(len(edges))
                                edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                                edges[str(len(edges)-1)].r_nei = e
                                del backup[i]
                                break
                            elif rr == 0:
                                end_edge = str(len(edges))
                                edges[e].r_nei = str(len(edges))
                                edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                                edges[str(len(edges)-1)].l_nei = e
                                del backup[i]
                                break
                
            

                boundary = edges[start_edge].get_loop(edges)

                # ax.plot(np.array(boundary)[:,0],np.array(boundary)[:,1],color='blue')
                # plt.show()

            if len(backup)!=0:
                for i in range(len(backup)-1,-1,-1):
                    edge = backup[i]
                    for e in [start_edge,end_edge]:
                        ll = np.linalg.norm(np.array(edge[0])-edges[e].l)
                        lr = np.linalg.norm(np.array(edge[0])-edges[e].r)
                        rl = np.linalg.norm(np.array(edge[1])-edges[e].l)
                        rr = np.linalg.norm(np.array(edge[1])-edges[e].r)
                        if ll == 0:
                            start_edge = str(len(edges))
                            edges[e].l_nei = str(len(edges))
                            edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                            edges[str(len(edges)-1)].r_nei = e
                            del backup[i]
                            break
                        elif lr == 0:
                            end_edge = str(len(edges))
                            edges[e].r_nei = str(len(edges))
                            edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                            edges[str(len(edges)-1)].l_nei = e
                            del backup[i]
                            break
                        elif rl == 0:
                            start_edge = str(len(edges))
                            edges[e].l_nei = str(len(edges))
                            edges[str(len(edges))] = Node(l=np.array(edge[0]),r=np.array(edge[1]))
                            edges[str(len(edges)-1)].r_nei = e
                            del backup[i]
                            break
                        elif rr == 0:
                            end_edge = str(len(edges))
                            edges[e].r_nei = str(len(edges))
                            edges[str(len(edges))] = Node(l=np.array(edge[1]),r=np.array(edge[0]))
                            edges[str(len(edges)-1)].l_nei = e
                            del backup[i]
                            break
            done.append(g)


            sequence = [-1,-1,-1,-1]
            values = [grids[g].r_nei,grids[g].t_nei,grids[g].l_nei,grids[g].b_nei]
            if grids[g].r_nei != None and grids[g].r_nei in done:
                sequence[0] = np.where(np.array(done)==grids[g].r_nei)[0][-1]
            elif grids[g].r_nei != None and grids[g].r_nei not in done:
                pass
            else:
                sequence[0] = np.Inf
            
            if grids[g].t_nei != None and grids[g].t_nei in done:
                sequence[1] = np.where(np.array(done)==grids[g].t_nei)[0][-1]
            elif grids[g].t_nei != None and grids[g].t_nei not in done:
                pass
            else:
                sequence[1] = np.Inf

            if grids[g].l_nei != None and grids[g].l_nei in done:
                sequence[2] = np.where(np.array(done)==grids[g].l_nei)[0][-1]
            elif grids[g].l_nei != None and grids[g].l_nei not in done:
                pass
            else:
                sequence[2] = np.Inf
            
            if grids[g].b_nei != None and grids[g].b_nei in done:
                sequence[3] = np.where(np.array(done)==grids[g].b_nei)[0][-1]
            elif grids[g].b_nei != None and grids[g].b_nei not in done:
                pass
            else:
                sequence[3] = np.Inf

            seq_ind = sequence.index(min(sequence))
            g1 = values[seq_ind]
            
            while g1 not in grids:
                sequence[sequence.index(min(sequence))] = np.Inf
                seq_ind = sequence.index(min(sequence))
                g1 = values[seq_ind]

            if g1 == grids[g].r_nei:
                g = grids[g].r_nei
            elif g1 == grids[g].t_nei:
                g = grids[g].t_nei
            elif g1 == grids[g].l_nei:
                g = grids[g].l_nei
            elif g1 == grids[g].b_nei:
                g = grids[g].b_nei

            if np.linalg.norm(edges[start_edge].l-edges[str(int(end_edge))].r)==0:
                break
        boundary = np.concatenate((boundary,np.array([boundary[0]])),axis=0)
        # ax.clear()
        # ax.plot(np.array(boundary)[:,0],np.array(boundary)[:,1],color='blue')
        # plt.show()
        # ax.clear()
        ax.plot(boundary[:,0],boundary[:,1],color = 'black')
        ax.scatter(boundary[:,0],boundary[:,1],color = 'black')
        plt.show()
        return boundary

    def Cut(p,c,ax):
        C_tr = p[str(c)].tr 
        
        boundary = polygon_boundary(p,ax)

        def get_neighbour(cel,loc,s):
            if s == 'r':
                return cel[loc].r_nei
            elif s == 'l':
                return cel[loc].l_nei
            elif s == 't':
                return cel[loc].t_nei
            elif s == 'b':
                return cel[loc].b_nei
        
        def set_neighbour(cells,loc,s,nei):
            if s == 'r':
                cells[loc].r_nei = nei
            elif s == 'l':
                cells[loc].l_nei = nei
            elif s == 't':
                cells[loc].t_nei = nei
            elif s == 'b':
                cells[loc].b_nei = nei
            return cells
        
        def get_area_vertices(array):
            array = np.array(array)
            done = 0
            pts = len(array)-2
            while not done:
                if pts<len(array)-1:
                    this = array[pts]
                    pre = array[pts-1]
                    post = array[pts+1]

                    if (pre[0] == this[0] == post[0]) or (pre[1] == this[1] == post[1]):
                        array = np.delete(array,pts,axis=0)
                pts -= 1
                if pts<1:
                    done = 1
            return array

        def get_opposite_side(s):
            if s == 'r':
                return 'l'
            elif s == 'l':
                return 'r'
            elif s == 't':
                return 'b'
            elif s == 'b':
                return 't'
        
        def get_perpendicular_sides(s):
            if s == 'r' or s == 'l':
                return ['t','b']
            elif s == 't' or s == 'b':
                return ['r','l']

        def get_region_vertex(i,j):
            if i+j == 'bl' or j+i == 'bl':
                return 'bl'
            elif i+j == 'br' or j+i == 'br':
                return 'br'
            elif i+j == 'tl' or j+i == 'tl':
                return 'tl'
            elif i+j == 'tr' or j+i == 'tr':
                return 'tr'
        
        def get_corner(i,j,grid):
            if i+j == 'bl' or j+i == 'bl':
                return grid.bl
            elif i+j == 'br' or j+i == 'br':
                return grid.br
            elif i+j == 'tl' or j+i == 'tl':
                return grid.tl
            elif i+j == 'tr' or j+i == 'tr':
                return grid.tr

        def check_presence_inside(i,j,vm,p2,C):
            if i+j == 'bl' or j+i == 'bl':
                a = p2[0]>vm[0] and p2[1]>vm[1] and p2[0]<C[0] and p2[1]<C[1]
                b = p2[0]>=vm[0] and p2[1]>vm[1] and p2[0]<C[0] and p2[1]<C[1]
                c = p2[0]>vm[0] and p2[1]>=vm[1] and p2[0]<C[0] and p2[1]<C[1]
                d = p2[0]>vm[0] and p2[1]>vm[1] and p2[0]<=C[0] and p2[1]<C[1]
                e = p2[0]>vm[0] and p2[1]>vm[1] and p2[0]<C[0] and p2[1]<=C[1]
                return a or b or c or d or e
            elif i+j == 'br' or j+i == 'br':
                a = p2[0]<vm[0] and p2[1]>vm[1] and p2[0]>C[0] and p2[1]<C[1]
                b = p2[0]<=vm[0] and p2[1]>vm[1] and p2[0]>C[0] and p2[1]<C[1]
                c = p2[0]<vm[0] and p2[1]>=vm[1] and p2[0]>C[0] and p2[1]<C[1]
                d = p2[0]<vm[0] and p2[1]>vm[1] and p2[0]>=C[0] and p2[1]<C[1]
                e = p2[0]<vm[0] and p2[1]>vm[1] and p2[0]>C[0] and p2[1]<=C[1]
                return a or b or c or d or e
            elif i+j == 'tl' or j+i == 'tl':
                a = p2[0]>vm[0] and p2[1]<vm[1] and p2[0]<C[0] and p2[1]>C[1]
                b = p2[0]>=vm[0] and p2[1]<vm[1] and p2[0]<C[0] and p2[1]>C[1]
                c = p2[0]>vm[0] and p2[1]<=vm[1] and p2[0]<C[0] and p2[1]>C[1]
                d = p2[0]>vm[0] and p2[1]<vm[1] and p2[0]<=C[0] and p2[1]>C[1]
                e = p2[0]>vm[0] and p2[1]<vm[1] and p2[0]<C[0] and p2[1]>=C[1]
                return a or b or c or d or e
            elif i+j == 'tr' or j+i == 'tr':
                a = p2[0]<vm[0] and p2[1]<vm[1] and p2[0]>C[0] and p2[1]>C[1]
                b = p2[0]<=vm[0] and p2[1]<vm[1] and p2[0]>C[0] and p2[1]>C[1]
                c = p2[0]<vm[0] and p2[1]<=vm[1] and p2[0]>C[0] and p2[1]>C[1]
                d = p2[0]<vm[0] and p2[1]<vm[1] and p2[0]>=C[0] and p2[1]>C[1]
                e = p2[0]<vm[0] and p2[1]<vm[1] and p2[0]>C[0] and p2[1]>=C[1]
                return a or b or c or d or e

        def get_v_m_tilda(i,j,p,start,stop,start_vertex):
            grid_list = []

            if isinstance(stop, type(None)):
                next_ = get_neighbour(p,start,i)
                grid_list.append(start)
                x = p[grid_list[-1]].get_x()
                y = p[grid_list[-1]].get_y()
                ax.plot(x,y)
                while next_ != stop:
                    grid_list.append(next_)
                    next_ = get_neighbour(p,next_,i)
                    x = p[grid_list[-1]].get_x()
                    y = p[grid_list[-1]].get_y()
                    ax.plot(x,y)
                
                next_ = get_neighbour(p,grid_list[-1],j)
                while next_ != stop:
                    grid_list.append(next_)
                    next_ = get_neighbour(p,next_,j)
                    x = p[grid_list[-1]].get_x()
                    y = p[grid_list[-1]].get_y()
                    ax.plot(x,y)
            else:
                next_nei = get_neighbour(p,start,i)
                grid_list.append(start)
                x = p[grid_list[-1]].get_x()
                y = p[grid_list[-1]].get_y()
                ax.plot(x,y)
                if next_nei != None:
                    while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                        grid_list.append(next_nei)
                        next_nei = get_neighbour(p,next_nei,i)
                        x = p[grid_list[-1]].get_x()
                        y = p[grid_list[-1]].get_y()
                        ax.plot(x,y)
                        if next_nei == None:
                            break
            
                next_nei = get_neighbour(p,grid_list[-1],j)
                if next_nei != None:
                    while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                        grid_list.append(next_nei)
                        next_nei = get_neighbour(p,next_nei,j)
                        x = p[grid_list[-1]].get_x()
                        y = p[grid_list[-1]].get_y()
                        ax.plot(x,y)
                        if next_nei == None:
                            break

            v_m_tilda = get_corner(i,j,p[grid_list[-1]])
            stop = v_m_tilda
            side = get_opposite_side(i)
            next_nei = get_neighbour(p,grid_list[-1],side)
            if next_nei != None:
                while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                    grid_list.append(next_nei)
                    x = p[grid_list[-1]].get_x()
                    y = p[grid_list[-1]].get_y()
                    ax.plot(x,y)
                    side_in = get_opposite_side(j)
                    next_nei_in = get_neighbour(p,grid_list[-1],side_in)
                    if next_nei_in != None:
                        while check_presence_inside(i,j,stop,p[next_nei_in].centroid,start_vertex):
                            grid_list.append(next_nei_in)
                            next_nei_in = get_neighbour(p,next_nei_in,side_in)
                            x = p[grid_list[-1]].get_x()
                            y = p[grid_list[-1]].get_y()
                            ax.plot(x,y)
                            if next_nei_in == None:
                                break
                    next_nei = get_neighbour(p,next_nei,side)
                    if next_nei == None:
                        break

            return v_m_tilda, grid_list

        region = get_area_vertices(boundary)
        sides = ['r','l','t','b']
        region_vertex = {'bl': [], 'br': [], 'tl': [], 'tr': []}
        for r in region:
            if r[0]>C_tr[0] and r[1]>C_tr[1]:
                region_vertex['tr'].append(r)
            elif r[0]<C_tr[0] and r[1]>C_tr[1]:
                region_vertex['tl'].append(r)
            elif r[0]>C_tr[0] and r[1]<C_tr[1]:
                region_vertex['br'].append(r)
            elif r[0]<C_tr[0] and r[1]<C_tr[1]:
                region_vertex['bl'].append(r)
        cut_success = {}
        # ax.scatter(C_tr[0],C_tr[1])
        # plt.show()
        for i in sides:
            both = get_perpendicular_sides(i)
            for j in both:
                key = i+j
                C_tr_grid = str(c)
                if key == 'lt' or key == 'tl':
                    C_tr_grid = p[C_tr_grid].t_nei
                elif key == 'rt' or key == 'tr':
                    C_tr_grid = p[p[C_tr_grid].t_nei].r_nei
                elif key == 'rb' or key == 'br':
                    C_tr_grid = p[C_tr_grid].r_nei

                
                vertex_inside = True
                v_m_tilda, grid_list = get_v_m_tilda(i,j,p,C_tr_grid,None,C_tr)
                stop_grid = v_m_tilda
                region_vertices = region_vertex[get_region_vertex(i,j)]
                count = -1
                while vertex_inside:
                    vertex_inside_list = []
                    for r in region_vertices:
                        vertex_inside_list.append(check_presence_inside(i,j,v_m_tilda,r,C_tr))

                    vertex_inside = np.sum(vertex_inside_list)
                    stop_grid = v_m_tilda
                    v_m_tilda, grid_list = get_v_m_tilda(i,j,p,C_tr_grid,stop_grid,C_tr)

                    if len(grid_list) == 1:    
                        marker = [0,0,0,0]  # ['r','t','l','b']
                        if p[grid_list[-1]].l_nei is not None:
                            marker[2] = 1
                        if p[grid_list[-1]].r_nei is not None:
                            marker[0] = 1
                        if p[grid_list[-1]].t_nei is not None:
                            marker[1] = 1
                        if p[grid_list[-1]].b_nei is not None:
                            marker[3] = 1                        
                    if not vertex_inside:
                        if len(grid_list) == 1:
                            if (i+j == 'bl' or j+i =='bl') and (marker[3]+marker[2] == 2):
                                grid_list = []
                            elif (i+j == 'br' or j+i =='br') and (marker[0]+marker[3] == 2):
                                grid_list = []
                            elif (i+j == 'tl' or j+i =='tl') and (marker[1]+marker[2] == 2):
                                grid_list = []
                            elif (i+j == 'tr' or j+i =='tr') and (marker[1]+marker[0] == 2):
                                grid_list = []
                        break
                    elif vertex_inside:
                        count += 1
                        v_m_tilda = region_vertices[count]

                vertex_exists = False
                corner_points = np.array([C_tr,np.array([C_tr[0],v_m_tilda[1]]),np.array([v_m_tilda[0],C_tr[1]])])
                for r in region:
                    for point in corner_points:
                        if np.linalg.norm(r-point)==0:
                            vertex_exists = True
                            break
                    if vertex_exists:
                        break
                if not vertex_exists:
                    cut_success[key] = np.unique(grid_list)
                else:
                    cut_success[key] = []
        
        cutting = False
        cutting_loc = None
        while len(cut_success) != 0 and cutting != True:
            checker = False
            random_ = np.random.randint(0,len(cut_success))
            for i,k in enumerate(cut_success):
                if i == random_ and len(cut_success[k])>0:
                    cutting = True
                    cutting_loc = k
                    break
                elif i == random_ and len(cut_success[k])==0:
                    checker = True
                    break
            if checker:
                del cut_success[k]
        
        if cutting_loc is not None:
            for g in cut_success[cutting_loc]:
                if g in cut_success[cutting_loc]:
                    if p[g].t_nei != None:
                        p[p[g].t_nei].b_nei = None
                    if p[g].b_nei != None:
                        p[p[g].b_nei].t_nei = None
                    if p[g].r_nei != None:
                        p[p[g].r_nei].l_nei = None
                    if p[g].l_nei != None:
                        p[p[g].l_nei].r_nei = None

                    # x = p[g].get_x()
                    # y = p[g].get_y()
                    # ax.plot(x,y,color='red')
                    # plt.show()
                    del p[g]

        return p,cutting

    def Inflate(p,c,ax):
        C_tr = p[str(c)].tr
        
        # ax.clear()            
        # for pp in p:
        #     x = p[pp].get_x()
        #     y = p[pp].get_y()
        #     ax.plot(x,y,color='red')
        #     plt.show()
        # x = p[str(c)].get_x()
        # y = p[str(c)].get_y()
        # ax.plot(x,y,color='green')
        # plt.show()

        new_cells = {}

        for ind in range(len(p)):
            grid = list(p.items())[ind][0]
            #   Quadrant marking
            # x = p[str(grid)].get_x()
            # y = p[str(grid)].get_y()
            # ax.plot(x,y,color='black')
            # plt.show()
            marker = np.zeros((4,4))
            grid_corners = p[grid].get_corners()
            for vertex in range(len(grid_corners)):
                if grid_corners[vertex,0] < C_tr[0] and grid_corners[vertex,1] < C_tr[1]:
                    marker[vertex,0] = 1
                elif grid_corners[vertex,0] < C_tr[0] and grid_corners[vertex,1] > C_tr[1]:
                    marker[vertex,3] = 1
                elif grid_corners[vertex,0] > C_tr[0] and grid_corners[vertex,1] < C_tr[1]:
                    marker[vertex,1] = 1
                elif grid_corners[vertex,0] > C_tr[0] and grid_corners[vertex,1] > C_tr[1]:
                    marker[vertex,2] = 1
                else:
                    pass

            #   Shifting  and Inserting cells
            a = None

            if np.sum(marker[:,3]) == 4:
                # NW
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,1] += 10

            if np.sum(marker[:,1]) == 4:
                # SE
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,0] += 10
            
            if np.sum(marker[:,2]) >= 1:
                # NE
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex] = grid_corners[vertex] + np.array([10,10])
            
            if grid_corners[2,1] == C_tr[1] and grid_corners[3,1] == C_tr[1] and grid_corners[0,0] >= C_tr[0]:
                #   +x bottom all
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,0] += 10

                t = Grid(bl = grid_corners[3], br = grid_corners[2], tl = grid_corners[3] + np.array([0,10]), tr = grid_corners[2] + np.array([0,10]))
                t.b_nei = grid
                if p[grid].t_nei is not None:
                    t.t_nei = p[grid].t_nei
                    p[t.t_nei].b_nei = str(int(list(p.items())[-1][0])+1)
                p[grid].t_nei = str(int(list(p.items())[-1][0])+1)
                new_cells[str(int(list(p.items())[-1][0])+1)] = t
                p.update(new_cells)
                # x = t.get_x()
                # y = t.get_y()
                # ax.plot(x,y,color='blue')
                # plt.show()


            if grid_corners[2,1] == C_tr[1] and grid_corners[3,1] == C_tr[1] and grid_corners[0,0] < C_tr[0]:
                #   -x bottom all

                a = Grid(bl = grid_corners[3], br = grid_corners[2], tl = grid_corners[3] + np.array([0,10]), tr = grid_corners[2] + np.array([0,10]))

                a.b_nei = grid
                if p[grid].t_nei is not None:
                    a.t_nei = p[grid].t_nei
                    p[p[grid].t_nei].b_nei = str(int(list(p.items())[-1][0])+1)
                p[grid].t_nei = str(int(list(p.items())[-1][0])+1)
                new_cells[str(int(list(p.items())[-1][0])+1)] = a
                p.update(new_cells)
                # x = a.get_x()
                # y = a.get_y()
                # ax.plot(x,y,color='blue')
                # plt.show()

            if grid_corners[1,0] == C_tr[0] and grid_corners[2,0] == C_tr[0] and grid_corners[3,1] <= C_tr[1]:
                #   -y left all

                a = Grid(bl = grid_corners[1], br = grid_corners[1] + np.array([10,0]), tl = grid_corners[2], tr = grid_corners[2] + np.array([10,0]))

                a.l_nei = grid
                if p[grid].r_nei is not None:
                    a.r_nei = p[grid].r_nei
                    p[p[grid].r_nei].l_nei = str(int(list(p.items())[-1][0])+1)
                p[grid].r_nei = str(int(list(p.items())[-1][0])+1)
                new_cells[str(int(list(p.items())[-1][0])+1)] = a
                p.update(new_cells)
                # x = a.get_x()
                # y = a.get_y()
                # ax.plot(x,y,color='blue')
                # plt.show()

                if grid_corners[3,1] == C_tr[1]:
                    a = Grid(bl = grid_corners[2], br = grid_corners[2] + np.array([10,0]), tl = grid_corners[2] + np.array([0,10]), tr = grid_corners[2] + np.array([10,10]))

                    a.l_nei = p[grid].t_nei
                    a.b_nei = p[grid].r_nei
                    p[p[grid].t_nei].r_nei = str(int(list(p.items())[-1][0])+1)
                    p[p[grid].r_nei].t_nei = str(int(list(p.items())[-1][0])+1)
                    
                    new_cells[str(int(list(p.items())[-1][0])+1)] = a
                    p.update(new_cells)
                    # x = a.get_x()
                    # y = a.get_y()
                    # ax.plot(x,y,color='blue')
                    # plt.show()

            if grid_corners[0,0] == C_tr[0] and grid_corners[3,0] == C_tr[0] and grid_corners[3,1] < C_tr[1]:
                #   -y right except first
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,0] += 10
            
            if grid_corners[0,1] == C_tr[1] and grid_corners[1,1] == C_tr[1] and grid_corners[2,0] < C_tr[0]:
                #   -x top except first
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,1] += 10

            if grid_corners[1,0] == C_tr[0] and grid_corners[2,0] == C_tr[0] and grid_corners[1,1] >= C_tr[1]:
                #   +y left excluding first cell
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,1] += 10

                a = Grid(bl = grid_corners[1], br = grid_corners[1] + np.array([10,0]), tl = grid_corners[2], tr = grid_corners[2] + np.array([10,0]))
                
                a.l_nei = grid
                if p[grid].r_nei is not None:
                    a.r_nei = p[grid].r_nei
                    p[p[grid].r_nei].l_nei = str(int(list(p.items())[-1][0])+1)
                p[grid].r_nei = str(int(list(p.items())[-1][0])+1)
                new_cells[str(int(list(p.items())[-1][0])+1)] = a
                p.update(new_cells)
                # x = a.get_x()
                # y = a.get_y()
                # ax.plot(x,y,color='blue')
                # plt.show()
            
            p[grid].bl = grid_corners[0]
            p[grid].br = grid_corners[1]
            p[grid].tr = grid_corners[2]
            p[grid].tl = grid_corners[3]
            p[grid].get_centroid()
        
            # ax.clear()            
            # for pp in p:
            #     x = p[pp].get_x()
            #     y = p[pp].get_y()
            #     ax.plot(x,y,color='red')
            #     plt.show()
            # x = p[str(c)].get_x()
            # y = p[str(c)].get_y()
            # ax.plot(x,y,color='green')
            # plt.show()
            # for pp in new_cells:
            #     x = new_cells[pp].get_x()
            #     y = new_cells[pp].get_y()
            #     ax.plot(x,y,color='blue')
            #     plt.show()
        
        p.update(new_cells)
        # ax.clear()            
        # for pp in p:
        #     x = p[pp].get_x()
        #     y = p[pp].get_y()
        #     ax.plot(x,y,color='red')
        #     plt.show()
        # x = p[str(c)].get_x()
        # y = p[str(c)].get_y()
        # ax.plot(x,y,color='green')
        # plt.show()
        # for pp in new_cells:
        #     x = new_cells[pp].get_x()
        #     y = new_cells[pp].get_y()
        #     ax.plot(x,y,color='blue')
        #     plt.show()

        for _,grid in enumerate(new_cells):
            # x = p[grid].get_x()
            # y = p[grid].get_y()
            # ax.plot(x,y,color='pink')
            # plt.show()
            if p[grid].l_nei is None:
                if p[grid].b_nei is not None:
                    if p[p[grid].b_nei].l_nei is not None:  
                        if p[p[p[grid].b_nei].l_nei].t_nei is not None:
                            p[grid].l_nei = p[p[p[grid].b_nei].l_nei].t_nei
                            if p[p[p[grid].b_nei].l_nei].t_nei is not None:
                                p[p[p[p[grid].b_nei].l_nei].t_nei].r_nei = grid
                        elif p[grid].t_nei is not None:
                            if p[p[grid].t_nei].l_nei is not None:
                                p[grid].l_nei = p[p[p[grid].t_nei].l_nei].b_nei
                                if p[p[p[grid].t_nei].l_nei].b_nei is not None:
                                    p[p[p[p[grid].t_nei].l_nei].b_nei].r_nei = grid
                elif p[grid].t_nei is not None:
                    if p[p[grid].t_nei].l_nei is not None:
                        p[grid].l_nei = p[p[p[grid].t_nei].l_nei].b_nei
                        if p[p[p[grid].t_nei].l_nei].b_nei is not None:
                            p[p[p[p[grid].t_nei].l_nei].b_nei].r_nei = grid
            if p[grid].r_nei is None:
                if p[grid].b_nei is not None:
                    if p[p[grid].b_nei].r_nei is not None:
                        if p[p[p[grid].b_nei].r_nei].t_nei is not None:
                            p[grid].r_nei = p[p[p[grid].b_nei].r_nei].t_nei
                            if p[p[p[grid].b_nei].r_nei].t_nei is not None:
                                p[p[p[p[grid].b_nei].r_nei].t_nei].l_nei = grid
                        elif p[grid].t_nei is not None:
                            if p[p[grid].t_nei].r_nei is not None:
                                p[grid].r_nei =  p[p[p[grid].t_nei].r_nei].b_nei
                                if p[p[p[grid].t_nei].r_nei].b_nei is not None:
                                    p[p[p[p[grid].t_nei].r_nei].b_nei].l_nei = grid
                elif p[grid].t_nei is not None:
                    if p[p[grid].t_nei].r_nei is not None:
                        p[grid].r_nei =  p[p[p[grid].t_nei].r_nei].b_nei
                        if p[p[p[grid].t_nei].r_nei].b_nei is not None:
                            p[p[p[p[grid].t_nei].r_nei].b_nei].l_nei = grid
            if p[grid].t_nei is None:
                if p[grid].r_nei is not None:
                    if p[p[grid].r_nei].t_nei is not None:
                        if p[p[p[grid].r_nei].t_nei].l_nei is not None:
                            p[grid].t_nei = p[p[p[grid].r_nei].t_nei].l_nei
                            if p[p[p[grid].r_nei].t_nei].l_nei is not None:
                                p[p[p[p[grid].r_nei].t_nei].l_nei].b_nei = grid
                        elif p[grid].l_nei is not None:
                            if p[p[grid].l_nei].t_nei is not None:
                                p[grid].t_nei = p[p[p[grid].l_nei].t_nei].r_nei
                                if p[p[p[grid].l_nei].t_nei].r_nei is not None:
                                    p[p[p[p[grid].l_nei].t_nei].r_nei].b_nei = grid
                elif p[grid].l_nei is not None:
                    if p[p[grid].l_nei].t_nei is not None:
                        p[grid].t_nei = p[p[p[grid].l_nei].t_nei].r_nei
                        if p[p[p[grid].l_nei].t_nei].r_nei is not None:
                            p[p[p[p[grid].l_nei].t_nei].r_nei].b_nei = grid
            if p[grid].b_nei is None:
                if p[grid].r_nei is not None:
                    if p[p[grid].r_nei].b_nei is not None:
                        if p[p[p[grid].r_nei].b_nei].l_nei is not None:
                            p[grid].b_nei = p[p[p[grid].r_nei].b_nei].l_nei
                            if p[p[p[grid].r_nei].b_nei].l_nei is not None:
                                p[p[p[p[grid].r_nei].b_nei].l_nei].t_nei = grid
                        elif p[grid].l_nei is not None:
                            if p[p[grid].l_nei].b_nei is not None:
                                p[grid].b_nei = p[p[p[grid].l_nei].b_nei].r_nei
                                if p[p[p[grid].l_nei].b_nei].r_nei is not None:
                                    p[p[p[p[grid].l_nei].b_nei].r_nei].t_nei = grid
                elif p[grid].l_nei is not None:
                    if p[p[grid].l_nei].b_nei is not None:
                        p[grid].b_nei = p[p[p[grid].l_nei].b_nei].r_nei
                        if p[p[p[grid].l_nei].b_nei].r_nei is not None:
                            p[p[p[p[grid].l_nei].b_nei].r_nei].t_nei = grid
            # if p[grid].t_nei is not None:
            #     ax.plot([p[grid].centroid[0],p[p[grid].t_nei].centroid[0]],[p[grid].centroid[1],p[p[grid].t_nei].centroid[1]])
            # if p[grid].b_nei is not None:
            #     ax.plot([p[grid].centroid[0],p[p[grid].b_nei].centroid[0]],[p[grid].centroid[1],p[p[grid].b_nei].centroid[1]])
            # if p[grid].l_nei is not None:
            #     ax.plot([p[grid].centroid[0],p[p[grid].l_nei].centroid[0]],[p[grid].centroid[1],p[p[grid].l_nei].centroid[1]])
            # if p[grid].r_nei is not None:
            #     ax.plot([p[grid].centroid[0],p[p[grid].r_nei].centroid[0]],[p[grid].centroid[1],p[p[grid].r_nei].centroid[1]])
            # plt.show()
        return p

    r = (num_vertices/2) - 2

    P = {'0': Grid(bl = np.array([50,50]), br = np.array([60,50]), tr = np.array([60,60]), tl = np.array([50,60]))}   # Unit square

    plt.ion()
    fig,ax = plt.subplots()
    fig1,ax1 = plt.subplots()
    while r>0:
        cut_success = False
        # polygon_boundary(P,ax)
        while not cut_success:
            p_trial = cp.copy(P)
            random_c = np.random.randint(0,len(p_trial))
            random_c = int(list(p_trial.items())[random_c][0])
            p_trial = Inflate(p_trial,random_c,ax)
            p_trial, cut_success = Cut(p_trial,random_c,ax1)
        P = p_trial
        ax1.clear()
        polygon_boundary(P,ax1)
        r -= 1

Inflate_Cut_algorithm(20)
# U and Z method
class Decomposition2:
    
    def __init__(self):
        self.grid_pts = None
        self.area_vertices = None
    
    def get_decomposed_graph(self,grid_points):
        self.grid_pts = np.copy(grid_points)
        self.area_vertices = self.get_area_vertices(grid_points)
        rectangles = []
        for i in range(len(self.area_vertices)+3):
            if i == 20:
                print(i)
            pt_1 = self.area_vertices[int((i)%len(self.area_vertices))]
            pt_2 = self.area_vertices[int((i+1)%len(self.area_vertices))]
            pt_3 = self.area_vertices[int((i+2)%len(self.area_vertices))]
            pt_4 = self.area_vertices[int((i+3)%len(self.area_vertices))]
            vector_1 = np.array(pt_2) - np.array(pt_1)
            vector_2 = np.array(pt_3) - np.array(pt_4)
            vec_1_norm = np.linalg.norm(vector_1)
            vec_2_norm = np.linalg.norm(vector_2)
            if vec_1_norm != 0:
                vector_1 = vector_1/vec_1_norm
            if vec_2_norm != 0:
                vector_2 = vector_2/vec_2_norm
            
            parallelism = np.dot(vector_1,vector_2)
            
            check_line = np.array(pt_2) - np.array(pt_3)

            com_inside = np.array([np.sum(self.area_vertices[:,0]),np.sum(self.area_vertices[:,1])])/len(self.area_vertices)

            com_line = com_inside - np.array(pt_3)

            com_cross = np.cross(com_line,check_line)

            if pt_2[0] == pt_3 [0] and parallelism != -1:
                min_side = min(vec_1_norm,vec_2_norm)
                if min_side == vec_1_norm:
                    rec = np.array([pt_1,pt_2,pt_3,np.array([pt_1[0],pt_3[1]]),pt_1])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
                else:
                    rec = np.array([np.array([pt_4[0],pt_2[1]]),pt_2,pt_3,pt_4,np.array([pt_4[0],pt_2[1]])])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
            elif pt_2[1] == pt_3[1] and parallelism != -1:
                min_side = min(vec_1_norm,vec_2_norm)
                if min_side == vec_1_norm:
                    rec = np.array([pt_1,pt_2,pt_3,np.array([pt_3[0],pt_1[1]]),pt_1])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
                else:
                    rec = np.array([np.array([pt_2[0],pt_4[1]]),pt_2,pt_3,pt_4,np.array([pt_2[0],pt_4[1]])])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)

        return rectangles

    def get_area_vertices(self,array):
        array = np.array(array)
        done = 0
        pts = len(array)-2
        while not done:
            if pts<len(array)-1:
                this = array[pts]
                pre = array[pts-1]
                post = array[pts+1]

                if (pre[0] == this[0] == post[0]) or (pre[1] == this[1] == post[1]):
                    array = np.delete(array,pts,axis=0)
            pts -= 1
            if pts<1:
                done = 1
        return array

# Modified U and Z method using 5 vertex
class Decomposition1:
    
    def __init__(self):
        self.grid_pts = None
        self.area_vertices = None
    
    def get_decomposed_graph(self,grid_points):
        self.grid_pts = np.copy(grid_points)
        self.area_vertices = self.get_area_vertices(grid_points)
        rectangles = []
        for i in range(len(self.area_vertices)+4):
            if i == 20:
                print(i)
            pt_1 = self.area_vertices[int((i)%len(self.area_vertices))]
            pt_2 = self.area_vertices[int((i+1)%len(self.area_vertices))]
            pt_3 = self.area_vertices[int((i+2)%len(self.area_vertices))]
            pt_4 = self.area_vertices[int((i+3)%len(self.area_vertices))]
            pt_5 = self.area_vertices[int((i+4)%len(self.area_vertices))]
            vector_1 = np.array(pt_2) - np.array(pt_1)
            vector_2 = np.array(pt_3) - np.array(pt_4)
            vec_1_norm = np.linalg.norm(vector_1)
            vec_2_norm = np.linalg.norm(vector_2)
            if vec_1_norm != 0:
                vector_1 = vector_1/vec_1_norm
            if vec_2_norm != 0:
                vector_2 = vector_2/vec_2_norm
            
            parallelism = np.dot(vector_1,vector_2)
            
            check_line = np.array(pt_2) - np.array(pt_3)

            com_inside = np.array([np.sum(self.area_vertices[:,0]),np.sum(self.area_vertices[:,1])])/len(self.area_vertices)

            com_line = com_inside - np.array(pt_3)

            com_cross = np.cross(com_line,check_line)

            if pt_2[0] == pt_3 [0] and parallelism != -1:
                min_side = min(vec_1_norm,vec_2_norm)
                if min_side == vec_1_norm:
                    rec = np.array([pt_1,pt_2,pt_3,np.array([pt_1[0],pt_3[1]]),pt_1])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
                else:
                    rec = np.array([np.array([pt_4[0],pt_2[1]]),pt_2,pt_3,pt_4,np.array([pt_4[0],pt_2[1]])])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
            elif pt_2[1] == pt_3[1] and parallelism != -1:
                min_side = min(vec_1_norm,vec_2_norm)
                if min_side == vec_1_norm:
                    rec = np.array([pt_1,pt_2,pt_3,np.array([pt_3[0],pt_1[1]]),pt_1])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)
                else:
                    rec = np.array([np.array([pt_2[0],pt_4[1]]),pt_2,pt_3,pt_4,np.array([pt_2[0],pt_4[1]])])
                    new_rec_com = np.array([np.sum(rec[:,0]),np.sum(rec[:,1])])/len(rec)
                    new_rec_com_line = new_rec_com - np.array(pt_3)
                    new_rec_com_cross = np.cross(new_rec_com_line,check_line)
                    check_parallism = np.dot(new_rec_com_cross,com_cross)
                    if check_parallism>=0:
                        rectangles.append(rec)

        return rectangles

    def get_area_vertices(self,array):
        array = np.array(array)
        done = 0
        pts = len(array)-2
        while not done:
            if pts<len(array)-1:
                this = array[pts]
                pre = array[pts-1]
                post = array[pts+1]

                if (pre[0] == this[0] == post[0]) or (pre[1] == this[1] == post[1]):
                    array = np.delete(array,pts,axis=0)
            pts -= 1
            if pts<1:
                done = 1
        return array

# Decomposition method
class Decomposition:
    
    def __init__(self):
        self.grid_pts = None
        self.area_vertices = None
    
    def get_decomposed_graph(self,grid_points):
        
        self.grid_pts = np.copy(grid_points)
        self.area_vertices = self.get_area_vertices(grid_points)

        H,V = self.get_bipartite(self.area_vertices,self.grid_pts)

        #H,V = self.get_absolute_bipartite(H,V)
        
        remaining_edges = self.get_remaining_edge(self.area_vertices,self.grid_pts,H,V)
        # print(H,V,remaining_edges)

        return H, V, remaining_edges

    def get_area_vertices(self,array):
        array = np.array(array)
        done = 0
        pts = len(array)-2
        while not done:
            if pts<len(array)-1:
                this = array[pts]
                pre = array[pts-1]
                post = array[pts+1]

                if (pre[0] == this[0] == post[0]) or (pre[1] == this[1] == post[1]):
                    array = np.delete(array,pts,axis=0)
            pts -= 1
            if pts<1:
                done = 1
        return array
    
    def get_bipartite(self,coordinates,main_arr):
        H_of_G = [[0,0]]
        V_of_G = [[0,0]] 

        action = [0,0,0,0]

        angles = self.get_angle(coordinates)
        coordinates_less_1 = coordinates[:-1]
        coordinates_less_1_len = len(coordinates_less_1)

        for i,point in enumerate(coordinates[:-1]):
            x1 = coordinates[i+1][0]
            y1 = coordinates[i+1][1]

            x_dif = x1 - coordinates[i][0]
            y_dif = y1 - coordinates[i][1]
            arr = np.append(coordinates_less_1[i:],coordinates_less_1[:-(coordinates_less_1_len-(i))],axis=0)
            arr1 = arr[1:]
            arr2 = arr[:-1]
            angle_ = np.append(angles[i:],angles[:-(coordinates_less_1_len-(i))],axis=0)
            angle_1 = angle_[1:]
            angle_2 = angle_[:-1]

            if y_dif >= 1:
                action[up] = 1
                temp1 = self.get_edge_bipartite(arr1,up,angle_1,main_arr)
                temp2 = self.get_edge_bipartite(arr2,down,angle_2,main_arr)
                if temp1 != None:
                    V_of_G = np.append(V_of_G,temp1,axis=0)
                elif temp2 != None:
                    V_of_G = np.append(V_of_G,temp2,axis=0)
            elif y_dif <= -1:
                action[down] = 1
                temp1 = self.get_edge_bipartite(arr1,down,angle_1,main_arr)
                temp2 = self.get_edge_bipartite(arr2,up,angle_2,main_arr)
                if temp1 != None:
                    V_of_G = np.append(V_of_G,temp1,axis=0)
                elif temp2 != None:
                    V_of_G = np.append(V_of_G,temp2,axis=0)
            elif x_dif >= 1:
                action[right] = 1
                temp1 = self.get_edge_bipartite(arr1,right,angle_1,main_arr)
                temp2 = self.get_edge_bipartite(arr2,left,angle_2,main_arr)
                if temp1 != None:
                    H_of_G = np.append(H_of_G,temp1,axis=0)
                elif temp2 != None:
                    H_of_G = np.append(H_of_G,temp2,axis=0)
            elif x_dif <= -1:
                action[left] = 1
                temp1 = self.get_edge_bipartite(arr1,left,angle_1,main_arr)
                temp2 = self.get_edge_bipartite(arr2,right,angle_2,main_arr)
                if temp1 != None:
                    H_of_G = np.append(H_of_G,temp1,axis=0)
                elif temp2 != None:
                    H_of_G = np.append(H_of_G,temp2,axis=0)
        H,V = self.unique_rows(H_of_G[1:]),self.unique_rows(V_of_G[1:])

        H,V = self.get_absolute_bipartite(H,V)

        return H,V
    
    def get_edge_bipartite(self,array,action,angle,main_arr):
        x,y = array[0][0],array[0][1]
        if action == up:
            for i,elem in enumerate(array[1:]):
                if elem[0] == x and elem[1] > y:
                    if angle[i+1] == 1:
                        if not self.get_obstacle(np.array(array),elem[0],elem[1],up,main_arr):
                            x1 = x
                            y1 = y
                            x2 = elem[0]
                            y2 = elem[1]
                            return [[x1,y1],[x2,y2]]
        elif action == right:
            for i,elem in enumerate(array[1:]):
                if elem[1] == y and elem[0] > x:
                    if angle[i+1] == 1:
                        if not self.get_obstacle(np.array(array),elem[0],elem[1],right,main_arr):
                            x1 = x
                            y1 = y
                            x2 = elem[0]
                            y2 = elem[1]
                            return [[x1,y1],[x2,y2]]
        elif action == down:
            for i,elem in enumerate(array[1:]):
                if elem[0] == x and elem[1] < y:
                    if angle[i+1] == 1:
                        if not self.get_obstacle(np.array(array),elem[0],elem[1],down,main_arr):
                            x1 = x
                            y1 = y
                            x2 = elem[0]
                            y2 = elem[1]
                            return [[x1,y1],[x2,y2]]
        elif action == left:
            for i,elem in enumerate(array[1:]):
                if elem[1] == y and elem[0] < x:
                    if angle[i+1] == 1:
                        if not self.get_obstacle(np.array(array),elem[0],elem[1],left,main_arr):
                            x1 = x
                            y1 = y
                            x2 = elem[0]
                            y2 = elem[1]
                            return [[x1,y1],[x2,y2]]
        return None
    
    def get_obstacle(self,array,x2,y2,action,main_arr):
        obstacle = 0
        x = array[0][0]
        y = array[0][1]
        
        if action == up or action == down:
            low = min(y,y2)
            high = max(y,y2)
            for j in range(low+1,high):
                for i in main_arr[1:]:
                    if np.array_equal(np.array([x,j]),np.array(i)):
                        obstacle = 1                        
                        break

        if action == right or action == left:
            low = min(x,x2)
            high = max(x,x2)
            for j in range(low+1,high):
                for i in main_arr[1:]:
                    if np.array_equal(np.array([j,y]),np.array(i)):
                        obstacle = 1
                        break
        return obstacle
    
    def get_absolute_bipartite(self,H,V):
        V_return = [[0,0]]
        for i in range(0,len(V),2):
            high_V = max(V[i][1],V[i+1][1])
            low_V = min(V[i][1],V[i+1][1])
            V_x = V[i][0]
            flag = 0
            for j in range(0,len(H),2):
                high_H = max(H[j][0],H[j+1][0])
                low_H = min(H[j][0],H[j+1][0])
                H_y = H[j][1]

                if H_y < high_V and H_y > low_V:

                    if V_x < high_H and V_x > low_H:

                        flag = 1
            if flag == 0:
                V_return = (np.append(V_return,[V[i],V[i+1]],axis=0)).tolist()
        return H,list(V_return[1:])
    
    def get_remaining_edge(self,array,main_arr,H,V):
        main_H = [[0,0]]
        main_V = [[0,0]]
        for i in range(0,len(H),2):
            low = min(H[i][0],H[i+1][0])
            high = max(H[i][0],H[i+1][0])
            
            for j in range(low,high+1):
                main_H = list(np.append(main_H,[[j,H[i][1]]],axis=0))
        main_H = main_H[1:]
        for i in range(0,len(V),2):
            low = min(V[i][1],V[i+1][1])
            high = max(V[i][1],V[i+1][1])
            
            for j in range(low,high+1):
                main_V = list(np.append(main_V,[[V[i][0],j]],axis=0))
        main_V = main_V[1:]
        
        all_points = main_arr
        if len(main_H) > 0:
            all_points = list(np.append(all_points,main_H,axis=0))
        if len(main_V) > 0:
            all_points = list(np.append(all_points,main_V,axis=0))
        
        all_points = np.unique(all_points,axis=0)
        
        angles = self.get_angle(array)
        
        remaining_point = [[0,0]]
        for i,elem in enumerate(array[:-1]):
            edge_flag = 0
            if angles[i] == 1:
                
                for point in H:
                    if np.array_equal(point,elem):
                        edge_flag = 1
                        break
                if edge_flag == 0:
                    
                    for point in V:
                        if np.array_equal(point,elem):
                            edge_flag = 1
                            break
                if edge_flag == 0:
                    remaining_point = list(np.append(remaining_point,[elem],axis=0))
                    
        remaining_point = remaining_point[1:]
        extra_edge = [[0,0]]

        for point in remaining_point:
#             print(remaining_point)
#             print(extra_edge)
            action = -1
            for i,elem in enumerate(main_arr):

                if np.array_equal(point,elem):
                    pre = main_arr[i-1]
                    post = main_arr[i+1]
                    if pre[1] < elem[1]:
                        action = up
                        break
                    elif pre[1] > elem[1]:
                        action = down
                        break
                    elif pre[1] == elem[1]:
                        if post[1] > elem[1]:
                            action = down
                        else:
                            action = up

            if action == up:
                extra_edge = np.append(extra_edge,[point],axis=0)
                yi = point[1]
                xi = point[0]
                exit = 0
                while(1):
                    yi+=1
                    #print(yi)
                    for j in all_points:
                        if np.array_equal(j,[xi,yi]):

                            extra_edge = np.append(extra_edge,[j],axis=0)
                            exit = 1
                            break
                        else:
                            pass
                    if exit == 1:
                        break

            if action == down:
                extra_edge = np.append(extra_edge,[point],axis=0)
                yi = point[1]
                xi = point[0]
                exit = 0
                while(1):
                    yi-=1
                    #print(xi,yi)
                    for j in all_points:
                        if np.array_equal(j,[xi,yi]):

                            extra_edge = np.append(extra_edge,[j],axis=0)
                            exit = 1
                            break
                        else:
                            pass
                    if exit == 1:
                        break

        extra_edge = extra_edge[1:]
        return list(extra_edge)

    def unique_rows(self,a):
        #print(a)
        for i in range(len(a)):
            for j in range(i):
                if np.array_equal(a[j],a[i]):

                    a[i]=[-1,-1]
                    #print(a)
                    break
        return_a = [elem.tolist() for i,elem in enumerate(a) if not np.array_equal([-1,-1],elem)]
        #print('final_a',return_a)
        return return_a

    def get_angle(self,array):
        array = np.array(array)
        angle = np.zeros(len(array)-1)
        action = np.array([0,0])

        for i in range(len(array)-1):
            pre = array[i-1]
            this = array[i]
            post = array[i+1] 

            pre_diff = (this - pre)
            if pre_diff[0] >= 1:
                action[0] = right
            elif pre_diff[0] <= -1:
                action[0] = left
            elif pre_diff[1] >= 1:
                action[0] = up
            elif pre_diff[1] <= -1:
                action[0] = down

            post_diff = (post - this)
            if post_diff[0] >= 1:
                action[1] = right
            elif post_diff[0] <= -1:
                action[1] = left
            elif post_diff[1] >= 1:
                action[1] = up
            elif post_diff[1] <= -1:
                action[1] = down

            if (action[0]==right and action[1]==down) or (action[0]==left and action[1]==up) or (action[0]==down and action[1]==left) or (action[0]==up and action[1]==right):
                angle[i-1] = 1    # Concave
            else:
                angle[i-1] = 0    # Convex
        return angle

def plot_main_figure(height,width,inner_grid_height,inner_grid_width):
    x_scalar,y_scalar = 0,0
    x,y = [0],[0]

    x_scalar,y_scalar = width,0
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)

    # print(x,y)
    x_scalar,y_scalar = x_scalar,y_scalar + height-inner_grid_height
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar + width,y_scalar
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar,y_scalar + height
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar - width,y_scalar
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar,y_scalar - height + inner_grid_height
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar-width ,y_scalar
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)

    x_scalar,y_scalar = x_scalar ,y_scalar-height
    x,y = np.append(x,x_scalar),np.append(y,y_scalar)
    # print(x,y)


    plt.plot(x,y,color = 'black')

    dual_graph_vertices_1, rectangular_hilbert_curves_1 = make_grid(x[0],y[0],x[6],y[6],inner_grid_height,inner_grid_width)
    dual_graph_vertices_2, rectangular_hilbert_curves_2 = make_grid(x[2],y[2],x[4],y[4],inner_grid_height,inner_grid_width)
    return dual_graph_vertices_1, dual_graph_vertices_2, rectangular_hilbert_curves_1, rectangular_hilbert_curves_2

def make_grid(x1,y1,x2,y2,inner_grid_height,inner_grid_width):
    #assume points convention (0,0),(1,0),(1,1),(0,1)

    x_center,y_center = [],[]

    width = x2-x1
    height = y2-y1


    hor_patch = [ x1 + i*inner_grid_width for i in range((width//inner_grid_width))]
    ver_patch = [ y1 + i*inner_grid_height for i in range((height//inner_grid_height))]

    print(hor_patch)
    print(ver_patch)

    for i,row in enumerate(hor_patch):
        for j,column in enumerate(ver_patch):

            x_scalar,y_scalar = row,column
            x_temp,y_temp =  row,column

            x_scalar,y_scalar = x_scalar+inner_grid_width,y_scalar
            x_temp,y_temp = np.append(x_temp,x_scalar),np.append(y_temp,y_scalar)

            x_scalar,y_scalar = x_scalar,y_scalar+inner_grid_height
            x_temp,y_temp = np.append(x_temp,x_scalar),np.append(y_temp,y_scalar)

            x_scalar,y_scalar = x_scalar-inner_grid_width,y_scalar
            x_temp,y_temp = np.append(x_temp,x_scalar),np.append(y_temp,y_scalar)

            x_scalar,y_scalar = x_scalar,y_scalar-inner_grid_height
            x_temp,y_temp = np.append(x_temp,x_scalar),np.append(y_temp,y_scalar)

            string = ("({},{})".format(row+(inner_grid_width//2),column+(inner_grid_height//2)))
            x_cen,y_cen = row+(inner_grid_width/2),column+(inner_grid_height/2)

            x_center,y_center = np.append(x_center,x_cen),np.append(y_center,y_cen)

            plt.plot(x_temp,y_temp,color='black')
        
    plt.scatter(x_center,y_center,c='green',s=10)

    rectangular_hilbert_curves = get_space_fill_graph(x_center[0],y_center[0],x_center[-1],y_center[-1],inner_grid_height,inner_grid_width)

    return [x_center,y_center], rectangular_hilbert_curves

#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-2-Clause
# Copyright (c) 2018 Jakub Červený

def gilbert2d(width, height):
    """
    Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
    2D rectangular grids. Generates discrete 2D coordinates to fill a rectangle
    of size (width x height).
    """

    if width >= height:
        yield from generate2d(0, 0, width, 0, 0, height)
    else:
        yield from generate2d(0, 0, 0, height, width, 0)

def sgn(x):
    return -1 if x < 0 else (1 if x > 0 else 0)

def generate2d(x, y, ax, ay, bx, by):

    w = abs(ax + ay)
    h = abs(bx + by)

    #print("this is w",w)

    (dax, day) = (sgn(ax), sgn(ay)) # unit major direction
    (dbx, dby) = (sgn(bx), sgn(by)) # unit orthogonal direction

    if h == 1:
        # trivial row fill
        for i in range(0, int(w)):
            yield(x, y)
            (x, y) = (x + dax, y + day)
        return

    if w == 1:
        # trivial column fill
        for i in range(0, int(h)):
            yield(x, y)
            (x, y) = (x + dbx, y + dby)
        return

    (ax2, ay2) = (ax//2, ay//2)
    (bx2, by2) = (bx//2, by//2)

    w2 = abs(ax2 + ay2)
    h2 = abs(bx2 + by2)

    if 2*w > 3*h:
        if (w2 % 2) and (w > 2):
            # prefer even steps
            (ax2, ay2) = (ax2 + dax, ay2 + day)

        # long case: split in two parts only
        yield from generate2d(x, y, ax2, ay2, bx, by)
        yield from generate2d(x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by)

    else:
        if (h2 % 2) and (h > 2):
            # prefer even steps
            (bx2, by2) = (bx2 + dbx, by2 + dby)

        # standard case: one step up, one long horizontal, one step down
        yield from generate2d(x, y, bx2, by2, ax2, ay2)
        yield from generate2d(x+bx2, y+by2, ax, ay, bx-bx2, by-by2)
        yield from generate2d(x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby),
                              -bx2, -by2, -(ax-ax2), -(ay-ay2))

def get_space_fill_graph(x1,y1,x2,y2,inner_grid_height,inner_grid_width):
    height = (abs(y2-y1)//inner_grid_height) + 1
    width = (abs(x2-x1)//inner_grid_width) + 1

    print("height : ",height)
    print("weight : ",width)

    x_arr = [ x*inner_grid_width for x, y in gilbert2d(width, height)]
    y_arr = [ y*inner_grid_height for x, y in gilbert2d(width, height)]

    x,y = x_arr + x1,y_arr + y1 

    plt.plot(x,y,color='red',linewidth = 2)

    return [x,y]



if __name__ == "__main__":
    up = 0
    right = 1
    down = 2
    left = 3

    graph = 0
    edge = 1

    # arena_w = 3
    # arena_h = 3
    # grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))

    # edges = adjacency_mat(grids)

    D = Decomposition()
    array = np.array([[0,0],[1,0],[2,0],[2,1],[3,1],[4,1],[4,0],[5,0],[5,1],[5,2],[5,3],[4,3],[4,4],[5,4],[5,5],[4,5],[3,5],[3,4],[2,4],[2,5],[1,5],[1,4],[0,4],[0,3],[1,3],[2,3],[3,3],[3,2],[2,2],[1,2],[1,1],[0,1],[0,0]])
    area_vertices = D.get_area_vertices(array)


    

    plt.ion()
    get_plot(array,graph)
    plot_show()

    H,V, remaining_edges = D.get_decomposed_graph(array)

    print(H,V, remaining_edges)

    get_plot(array,graph)
    get_plot(H,edge)
    get_plot(V,edge)
    get_plot(remaining_edges,edge)
    plot_show()

    rectangles = np.array(D.get_decomposed_graph(array))
    
    for r in range(len(rectangles)-1,-1,-1):
        check = []
        for i in range(len(rectangles[0])):
            if np.linalg.norm(rectangles[r,i]-rectangles[r-1,i])==0:
                check.append(1)
        if np.sum(check)>1:
            rectangles = np.delete(rectangles,r,axis=0)
            rectangles = np.delete(rectangles,r-1,axis=0)

    fig,ax = plt.subplots()
    for r in rectangles:
        ax.plot(r[:,0],r[:,1])
    
    plt.ioff()
    plt.show()

    # height = 120
    # width = 120

    # inner_grid_width = 15
    # inner_grid_height = 10

    # DG1, DG2, RHC1, RHC2 = plot_main_figure(height,width,inner_grid_height,inner_grid_width)
    
    # k_max = len(DG1[0]) + len(DG2[0])
    # time = []
    # optimal_k = []
    # for t in range(1,k_max+1):
    #     optimal_searcher = []
    #     for k in range(k_max,0,-1):
    #         if (len(DG1[0])+len(DG2[0]))/k - t == 0.0:
    #             optimal_searcher.append(k)

    #     if len(optimal_searcher)!=0:
    #         time.append(t)
    #         optimal_k.append(min(optimal_searcher))

            

    # plt.plot(time,optimal_k)
    # plt.show()









