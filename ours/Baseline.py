#   Authors: Swadhin Agrawal

import matplotlib.pyplot as plt
import numpy as np
import copy as cp
import pickle as pkl
import os
from shapely.geometry import LineString

# sed = 19778167 # 5 spike 19778167: 38lg
sed = 27819661 # 7 spike 27819661 38lg
np.random.seed(sed)

# Random simple rectilinear polygon generator

class Edge:
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
        
        self.sfc_path_pre = None
        self.sfc_path_post = None
        self.centroid = (self.tl+self.bl+self.tr+self.br)/4.0

        self.robot_home = False
        self.passage = 0
        self.passage_l = 0
        self.passage_r = 0
        self.passage_t = 0
        self.passage_b = 0

        self.rectangle_num = None

        self.chance_of_intruder = 1
        self.cost_of_traversal = -np.Inf

        self.heuristics = 0#np.random.uniform(0,0.1)
        self.heuristics_2 = 0
        self.distance = 0
        self.total_cost = 0

        self.robo_path_pre = None

        self.height = None
        self.width = None

        self.body = None
    
    def get_centroid(self):
        self.centroid = (self.tl+self.bl+self.tr+self.br)/4.0

    def get_x(self):
        return np.array([self.bl[0],self.br[0],self.tr[0],self.tl[0],self.bl[0]])

    def get_y(self):
        return np.array([self.bl[1],self.br[1],self.tr[1],self.tl[1],self.bl[1]])

    def get_corners(self):
        return np.array([self.bl,self.br,self.tr,self.tl])

def Inflate_Cut_algorithm(num_vertices,g_size=10,ax=None):
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
                            edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                        else:
                            for e in [start_edge,end_edge]:
                                ll = np.linalg.norm(np.array(edge[0])-edges[e].l)
                                lr = np.linalg.norm(np.array(edge[0])-edges[e].r)
                                rl = np.linalg.norm(np.array(edge[1])-edges[e].l)
                                rr = np.linalg.norm(np.array(edge[1])-edges[e].r)
                                if ll == 0:
                                    start_edge = str(len(edges))
                                    edges[e].l_nei = str(len(edges))
                                    edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
                                    edges[str(len(edges)-1)].r_nei = e
                                    break
                                elif lr == 0:
                                    end_edge = str(len(edges))
                                    edges[e].r_nei = str(len(edges))
                                    edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                                    edges[str(len(edges)-1)].l_nei = e
                                    break
                                elif rl == 0:
                                    start_edge = str(len(edges))
                                    edges[e].l_nei = str(len(edges))
                                    edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                                    edges[str(len(edges)-1)].r_nei = e
                                    break
                                elif rr == 0:
                                    end_edge = str(len(edges))
                                    edges[e].r_nei = str(len(edges))
                                    edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
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
                                edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
                                edges[str(len(edges)-1)].r_nei = e
                                del backup[i]
                                break
                            elif lr == 0:
                                end_edge = str(len(edges))
                                edges[e].r_nei = str(len(edges))
                                edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                                edges[str(len(edges)-1)].l_nei = e
                                del backup[i]
                                break
                            elif rl == 0:
                                start_edge = str(len(edges))
                                edges[e].l_nei = str(len(edges))
                                edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                                edges[str(len(edges)-1)].r_nei = e
                                del backup[i]
                                break
                            elif rr == 0:
                                end_edge = str(len(edges))
                                edges[e].r_nei = str(len(edges))
                                edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
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
                            edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
                            edges[str(len(edges)-1)].r_nei = e
                            del backup[i]
                            break
                        elif lr == 0:
                            end_edge = str(len(edges))
                            edges[e].r_nei = str(len(edges))
                            edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                            edges[str(len(edges)-1)].l_nei = e
                            del backup[i]
                            break
                        elif rl == 0:
                            start_edge = str(len(edges))
                            edges[e].l_nei = str(len(edges))
                            edges[str(len(edges))] = Edge(l=np.array(edge[0]),r=np.array(edge[1]))
                            edges[str(len(edges)-1)].r_nei = e
                            del backup[i]
                            break
                        elif rr == 0:
                            end_edge = str(len(edges))
                            edges[e].r_nei = str(len(edges))
                            edges[str(len(edges))] = Edge(l=np.array(edge[1]),r=np.array(edge[0]))
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
        if ax is not None:
            ax.plot(boundary[:,0],boundary[:,1],color = 'black')
            # ax.scatter(boundary[:,0],boundary[:,1],color = 'black')
            plt.show()
        return boundary, edges, start_edge, end_edge

    def Cut(p,c,ax=None):
        C_tr = p[str(c)].tr 
        
        boundary, _,_,_ = polygon_boundary(p,ax)

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
                # x = p[grid_list[-1]].get_x()
                # y = p[grid_list[-1]].get_y()
                # ax.plot(x,y)
                while next_ != stop:
                    grid_list.append(next_)
                    next_ = get_neighbour(p,next_,i)
                    # x = p[grid_list[-1]].get_x()
                    # y = p[grid_list[-1]].get_y()
                    # ax.plot(x,y)
                
                next_ = get_neighbour(p,grid_list[-1],j)
                while next_ != stop:
                    grid_list.append(next_)
                    next_ = get_neighbour(p,next_,j)
                    # x = p[grid_list[-1]].get_x()
                    # y = p[grid_list[-1]].get_y()
                    # ax.plot(x,y)
            else:
                next_nei = get_neighbour(p,start,i)
                grid_list.append(start)
                # x = p[grid_list[-1]].get_x()
                # y = p[grid_list[-1]].get_y()
                # ax.plot(x,y)
                if next_nei != None:
                    while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                        grid_list.append(next_nei)
                        next_nei = get_neighbour(p,next_nei,i)
                        # x = p[grid_list[-1]].get_x()
                        # y = p[grid_list[-1]].get_y()
                        # ax.plot(x,y)
                        if next_nei == None:
                            break
            
                next_nei = get_neighbour(p,grid_list[-1],j)
                if next_nei != None:
                    while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                        grid_list.append(next_nei)
                        next_nei = get_neighbour(p,next_nei,j)
                        # x = p[grid_list[-1]].get_x()
                        # y = p[grid_list[-1]].get_y()
                        # ax.plot(x,y)
                        if next_nei == None:
                            break

            v_m_tilda = get_corner(i,j,p[grid_list[-1]])
            stop = v_m_tilda
            side = get_opposite_side(i)
            next_nei = get_neighbour(p,grid_list[-1],side)
            if next_nei != None:
                while check_presence_inside(i,j,stop,p[next_nei].centroid,start_vertex):
                    grid_list.append(next_nei)
                    # x = p[grid_list[-1]].get_x()
                    # y = p[grid_list[-1]].get_y()
                    # ax.plot(x,y)
                    side_in = get_opposite_side(j)
                    next_nei_in = get_neighbour(p,grid_list[-1],side_in)
                    if next_nei_in != None:
                        while check_presence_inside(i,j,stop,p[next_nei_in].centroid,start_vertex):
                            grid_list.append(next_nei_in)
                            next_nei_in = get_neighbour(p,next_nei_in,side_in)
                            # x = p[grid_list[-1]].get_x()
                            # y = p[grid_list[-1]].get_y()
                            # ax.plot(x,y)
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

    def Inflate(p,c,g_size=10,ax=None):
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
                    grid_corners[vertex,1] += g_size

            if np.sum(marker[:,1]) == 4:
                # SE
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,0] += g_size
            
            if np.sum(marker[:,2]) >= 1:
                # NE
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex] = grid_corners[vertex] + np.array([g_size,g_size])
            
            if grid_corners[2,1] == C_tr[1] and grid_corners[3,1] == C_tr[1] and grid_corners[0,0] >= C_tr[0]:
                #   +x bottom all
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,0] += g_size

                t = Grid(bl = grid_corners[3], br = grid_corners[2], tl = grid_corners[3] + np.array([0,g_size]), tr = grid_corners[2] + np.array([0,g_size]))
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

                a = Grid(bl = grid_corners[3], br = grid_corners[2], tl = grid_corners[3] + np.array([0,g_size]), tr = grid_corners[2] + np.array([0,g_size]))

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

                a = Grid(bl = grid_corners[1], br = grid_corners[1] + np.array([g_size,0]), tl = grid_corners[2], tr = grid_corners[2] + np.array([g_size,0]))

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
                    a = Grid(bl = grid_corners[2], br = grid_corners[2] + np.array([g_size,0]), tl = grid_corners[2] + np.array([0,g_size]), tr = grid_corners[2] + np.array([g_size,g_size]))

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
                    grid_corners[vertex,0] += g_size
            
            if grid_corners[0,1] == C_tr[1] and grid_corners[1,1] == C_tr[1] and grid_corners[2,0] < C_tr[0]:
                #   -x top except first
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,1] += g_size

            if grid_corners[1,0] == C_tr[0] and grid_corners[2,0] == C_tr[0] and grid_corners[1,1] >= C_tr[1]:
                #   +y left excluding first cell
                for vertex in range(len(grid_corners)):
                    grid_corners[vertex,1] += g_size

                a = Grid(bl = grid_corners[1], br = grid_corners[1] + np.array([g_size,0]), tl = grid_corners[2], tr = grid_corners[2] + np.array([g_size,0]))
                
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

    P = {'0': Grid(bl = np.array([50,50]), br = np.array([50+g_size,50]), tr = np.array([50+g_size,50+g_size]), tl = np.array([50,50+g_size]))}   # Unit square

    plt.ion()
    # fig,ax = plt.subplots()
    if ax == None:
        fig1,ax1 = plt.subplots()
    else:
        ax1 = ax
    while r>0:
        cut_success = False
        # polygon_boundary(P,ax)
        while not cut_success:
            p_trial = cp.copy(P)
            random_c = np.random.randint(0,len(p_trial))
            random_c = int(list(p_trial.items())[random_c][0])
            p_trial = Inflate(p_trial,random_c,g_size=g_size)

            p_trial, cut_success = Cut(p_trial,random_c)
        P = p_trial
        ax1.clear()
        polygon_boundary(P,ax1)
        r -= 1
    boundary,edges,start_edge,end_edge = polygon_boundary(P,ax)
    return boundary,edges,start_edge,end_edge

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

# Grid graph construction:

def is_grid_in_poly(polygon,grid,min_x,min_y,max_x,max_y):
    c = grid.centroid
    l = np.array([c[0]-(max_x-min_x),c[1]])
    r = np.array([c[0]+(max_x-min_x),c[1]])
    t = np.array([c[0],c[1]+(max_y-min_y)])
    b = np.array([c[0],c[1]-(max_y-min_y)])
    cl = LineString([tuple(c),tuple(l)])
    cr = LineString([tuple(c),tuple(r)])
    ct = LineString([tuple(c),tuple(t)])
    cb = LineString([tuple(c),tuple(b)])
    B_l = []
    B_r = []
    B_t = []
    B_b = []
    for v in range(-1,len(polygon)-1):
        e = LineString([polygon[v],polygon[v+1]])

        B_l.append(e.intersects(cl))
        B_r.append(e.intersects(cr))
        B_t.append(e.intersects(ct))
        B_b.append(e.intersects(cb))
    if np.sum(B_l)%2==1 or np.sum(B_r)%2==1 or np.sum(B_t)%2==1 or np.sum(B_b)%2==1:
        return 1
    else:
        return 0

def Grid_graph(polygon,grid_size,ax):
    min_x = min(polygon[:,0])
    max_x = max(polygon[:,0])
    min_y = min(polygon[:,1])
    max_y = max(polygon[:,1])
    grids = {}
    
    for i in range(min_x,max_x,grid_size):
        for j in range(min_y,max_y,grid_size):
            tl = np.array([i, j+grid_size])
            tr = np.array([i+grid_size, j+grid_size])
            bl = np.array([i, j])
            br = np.array([i+grid_size, j])
            g = Grid(tl=tl,bl=bl,tr=tr,br=br)
            if is_grid_in_poly(polygon,g,min_x,min_y,max_x,max_y):
                if str(bl - np.array([grid_size,0])) in grids:
                    g.l_nei = str(bl - np.array([grid_size,0]))
                    grids[str(bl - np.array([grid_size,0]))].r_nei = str(bl)
                    # ax.plot([g.centroid[0],grids[str(bl - np.array([grid_size,0]))].centroid[0]],[g.centroid[1],grids[str(bl - np.array([grid_size,0]))].centroid[1]],color='maroon',linewidth=0.5)
                if str(br) in grids:
                    g.r_nei = str(br)
                    grids[str(br)].l_nei = str(bl)
                    # ax.plot([g.centroid[0],grids[str(br)].centroid[0]],[g.centroid[1],grids[str(br)].centroid[1]],color='maroon',linewidth=0.5)
                if str(tl) in grids:
                    g.t_nei = str(tl)
                    grids[str(tl)].b_nei = str(bl)
                    # ax.plot([g.centroid[0],grids[str(tl)].centroid[0]],[g.centroid[1],grids[str(tl)].centroid[1]],color='maroon',linewidth=0.5)
                if str(bl - np.array([0,grid_size])) in grids:
                    g.b_nei = str(bl - np.array([0,grid_size]))
                    grids[str(bl - np.array([0,grid_size]))].t_nei = str(bl)
                    # ax.plot([g.centroid[0],grids[str(bl - np.array([0,grid_size]))].centroid[0]],[g.centroid[1],grids[str(bl - np.array([0,grid_size]))].centroid[1]],color='maroon',linewidth=0.5)
                grids[str(bl)] = g
                # x = grids[str(bl)].get_x()
                # y = grids[str(bl)].get_y()
                # ax.plot(x,y,color='green',linewidth=0.5)
                
    return grids

#   Path planning for searcher robots              

def Astar(grids,source,destination,fig = None,ax = None):
    def get_neighbour(cel,s):
        if s == 'r':
            return cel.r_nei
        elif s == 'l':
            return cel.l_nei
        elif s == 't':
            return cel.t_nei
        elif s == 'b':
            return cel.b_nei

    temp_grids = cp.copy(grids)

    if ax != None:
        plot_mesh(fig,ax,grids,title = '',without_bar=1)

    temp_grids[source].total_cost = 0
    temp_grids[source].distance = 0
    temp_grids[destination].total_cost = 0
    
    if ax != None:
        temp_grids[source].body = ax.plot(temp_grids[source].get_x(),temp_grids[source].get_y(),color = 'white')
        temp_grids[destination].body = ax.plot(temp_grids[destination].get_x(),temp_grids[destination].get_y(),color = 'green')

    open_nodes = [source]
    closed_nodes = []
    while len(open_nodes) > 0:
        
        current_node = open_nodes[0]
        current_index = 0
        for i,item in enumerate(open_nodes):
            if temp_grids[item].total_cost < temp_grids[current_node].total_cost:
                current_node = item
                current_index = i

        # Pop current off open list, add to closed list
        open_nodes.pop(current_index)
        closed_nodes.append(current_node)
        # Found the goal
        if current_node == destination or temp_grids[current_node].l_nei == destination or temp_grids[current_node].r_nei == destination or temp_grids[current_node].t_nei == destination or temp_grids[current_node].b_nei == destination:
            temp_grids[destination].robo_path_pre = current_node
            path = []
            current = destination
            while current is not None:
                path.append(current)
                current = temp_grids[current].robo_path_pre
            
            # for p in path:
            #     temp_grids[p].body = ax.plot(temp_grids[p].get_x(),temp_grids[p].get_y(),color = 'yellow',linewidth = 0.5)
            
            return path[::-1] # Return reversed path

        # Generate child paths

        children = []
        sides = ['l','r','t','b']
        for s in sides:
            nei = get_neighbour(temp_grids[current_node],s)
            if nei != source:
                children.append(nei)

        # Loop through children
        for child in children:
            if child in closed_nodes:
                continue
            elif child != None:
                # Assign the f, g, and h values
                temp_grids[child].distance = temp_grids[current_node].distance + 1
                temp_grids[child].robo_path_pre = current_node
                temp_grids[child].total_cost = temp_grids[current_node].total_cost + temp_grids[child].distance#+ temp_grids[child].heuristics# + 0.1*temp_grids[child].heuristics_2

            # Add the child to the open list
            if child not in open_nodes+closed_nodes and child != None:
                open_nodes.append(child)
    # if current_node == destination or temp_grids[current_node].l_nei == destination or temp_grids[current_node].r_nei == destination or temp_grids[current_node].t_nei == destination or temp_grids[current_node].b_nei == destination:
    temp_grids[destination].robo_path_pre = current_node
    path = []
    current = destination
    while current is not None:
        path.append(current)
        current = temp_grids[current].robo_path_pre
    
    # for p in path:
    #     temp_grids[p].body = ax.plot(temp_grids[p].get_x(),temp_grids[p].get_y(),color = 'yellow',linewidth = 0.5)
    
    return path[::-1] # Return reversed path

    # return path

def plot_mesh(fig,ax,grid,title,grid_mesh,without_bar=0):
    values_x = []
    values_y = []
    for g in grid:
        values_x = values_x + list(grid[g].get_x())
        values_y = values_y + list(grid[g].get_y())

    x = np.arange(min(values_x),max(values_x),grid_width)
    y = np.arange(min(values_y),max(values_y),grid_height)
    z = []
    for i in x:
        verti = []
        for j in y:
            if str(np.array([i, j])) in grid:
                verti.append(grid[str(np.array([i, j]))].heuristics)
            else:
                verti.append(0)
        z.append(verti)
    z = np.array(z).T
    # set_array
    if grid_mesh == None:
        grid_mesh = ax.pcolormesh(x + grid_width/2,y + grid_height/2,z,shading='auto',cmap='gist_heat_r')
        grid_mesh.set_clim(0,1)
        if not without_bar:
            cbar = fig.colorbar(grid_mesh,orientation='vertical')
            cbar.set_ticks(np.arange(0,1.1,0.1))

        # ax.set_aspect('equal', 'box')
        # ax.set_xticks([])
        # ax.set_yticks([])
        ax.set_title(title)
        return grid_mesh
    else:
        grid_mesh.set_array(z)

#   Setting up robot object

class Robot:
    def __init__(self,id,r_type,grids):
        self.id = id
        self.type = r_type

        self.past_loc = None
        self.present_loc = np.random.choice(list(grids))
        self.next_loc = np.random.choice(list(grids))
        while self.next_loc == self.present_loc:
            self.next_loc = np.random.choice(list(grids))
        self.body = None

        self.path_history = []
        if self.type == 's':
            self.path_ahead = Astar(grids,self.present_loc,self.next_loc)

        self.path_progressor = 0

        self.trajectory = []
        self.plott = None

    def update_trajectory(self,loc):
        self.trajectory.pop(0)
        self.trajectory = self.trajectory + [loc]

    def loc(self,nodes,bias=0):
        if bias:
            probabilities = []
            for g in nodes:
                probabilities.append(nodes[g].chance_of_intruder)

        def get_neighbour(cel,s):
            if s == 'r':
                return cel.r_nei
            elif s == 'l':
                return cel.l_nei
            elif s == 't':
                return cel.t_nei
            elif s == 'b':
                return cel.b_nei
        
        if self.type == 's':
            return np.random.choice(list(nodes))
        else:
            sides = ['l','r','t','b']
            intruder_next = np.random.choice(sides)
            next_loc = get_neighbour(nodes[self.present_loc],intruder_next)
            while next_loc == None:
                sides.pop(sides.index(intruder_next))
                intruder_next = np.random.choice(sides)
                next_loc = get_neighbour(nodes[self.present_loc],intruder_next)
            return next_loc
    def update(self,grid):
        self.body.set_offsets([[grid[self.present_loc].centroid[0],grid[self.present_loc].centroid[1]]])

def intruder_(grids,grid_w,grid_h,ax):
    intruder = Robot(-1,'i',grids)    #   Spawn Intruder
    if ax != None:
        intruder.body = ax.scatter([grids[intruder.present_loc].centroid[0]],[grids[intruder.present_loc].centroid[1]],color='red',s=25)
        ax.add_artist(intruder.body)
    for g in grids:
        grids[g].robo_path_pre = None
    return intruder

def searchers(num_robots,grid_w,grid_h,grid,ax = None):
    robots = [] #   Spawn Robots
    for i in range(num_robots):
        robots.append(Robot(i,'s',grid))
        if ax != None:
            robots[i].body = ax.scatter([grid[robots[i].present_loc].centroid[0]],[grid[robots[i].present_loc].centroid[1]],color='green',s=25)
            ax.add_artist(robots[i].body)
        for g in grid:
            grid[g].robo_path_pre = None
    return robots

def update_prob_costs(grids,present_locs):
    # for g in grids: #   Update the Costs and Probabilities of the grids
    for l in present_locs:
        # if str(l)==g:
        # grids[g].probability -= 0#1/(arena_h*arena_w)
        # if grids[g].probability<0:
        #     grids[g].probability=0
        # elif grids[g].probability>1:
        #     grids[g].probability=1
        grids[str(l)].heuristics += 0.05#50/len(grids)
        # else:
        #     grids[g].probability += 0#1/(arena_h*arena_w*len(present_locs))
    # return grids

def dump_mesh(grid):
    values_x = []
    values_y = []
    for g in grid:
        values_x = values_x + list(grid[g].get_x())
        values_y = values_y + list(grid[g].get_y())

    x = np.arange(min(values_x),max(values_x),grid_width)
    y = np.arange(min(values_y),max(values_y),grid_height)
    z = []
    for i in x:
        verti = []
        for j in y:
            if str([i, j]) in grid:
                verti.append(grid[str([i,j])].heuristics)
            else:
                verti.append(0)
        z.append(verti)
    z = np.array(z).T
    return x,y,z

# Single robot
static_intruder = 0
if static_intruder:
    path_ = os.path.realpath(os.path.dirname(__file__))
    # fig,ax = plt.subplots()
    # plt.box(False)
    # ax.set_aspect('equal')
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    ax = None
    itera = np.random.randint(50)
    print('itera :',itera)
    
    tile_sixe = 10
    grid_height = 5
    grid_width = 5
    
    # Random simple connected rectilinear polygon
    boundary,edges,start_edge,end_edge = Inflate_Cut_algorithm(20,g_size=tile_sixe)#,ax=ax)
    boundary1 = get_area_vertices(np.concatenate((boundary[2:],boundary[:3]),axis=0))
    # ax.clear()
    # ax.plot(boundary1[:,0],boundary1[:,1],color='indigo')
    
    grid_graph = Grid_graph(boundary1,grid_size=grid_height,ax=ax)
    
    # plt.ioff()
    # plt.show()

    data_his = []
    # plt.ion()
    grid_mesh = None
    for num_robots in range(1,153,3):
        avg_time = []
        data_runs = []
        print(num_robots)
        for runs in range(100):
            
            for g in grid_graph:
                grid_graph[g].robo_path_pre = None
                grid_graph[g].heuristics = 0
                grid_graph[g].distance = 0
                grid_graph[g].total_cost = 0

            t = 0   #   Initialize time
            search_state =  0 #   0 = not found, 1 = found


            robots = searchers(num_robots,grid_width,grid_height,grid_graph,ax)
            intruder = intruder_(grid_graph,grid_width,grid_height,ax)
            
            #   Robots path
            for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                robots[i].next_loc = intruder.present_loc

                for g in grid_graph:
                    grid_graph[g].robo_path_pre = None
                
                if robots[i].present_loc != robots[i].next_loc:
                    path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)
                else:
                    path = [robots[i].next_loc]

                robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                robots[i].path_ahead = path
                robots[i].path_progressor = 0
                
            while not search_state: #   While Intruder is not found
                for r in robots:
                    if intruder.present_loc == r.present_loc:
                        # intruder.body.set_visible(False)
                        # for rb in robots:
                        #     rb.body.set_visible(False)
                        #     if rb.plott!= None:
                        #         rb.plott[0].set_visible(False)
                        search_state = 1
                        avg_time.append(t)
                present_locs = [r.present_loc for r in robots]
                update_prob_costs(grid_graph,present_locs)
                # if grid_mesh is None:
                #     grid_mesh = plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=0)
                # else:    
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=1)

                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_progressor = 0
                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    robots[i].path_progressor += 1

                    if len(robots[i].trajectory)<5:
                        robots[i].trajectory.append([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                        # if len(robots[i].trajectory)==4:
                        #     robots[i].plott = ax.plot(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                    else:
                        robots[i].update_trajectory([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     robots[i].plott[0].set_data(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1])
                    # robots[i].body.set_offsets([np.array(robots[i].trajectory)[-1,0],np.array(robots[i].trajectory)[-1,1]])

                # Intruder update
                intruder.next_loc = intruder.present_loc
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc

                t += 1
                # if runs==1 and t==620:
                #     path = os.getcwd()
                #     plt.savefig(path+'/results/1RIS_s_bench.pdf',format = "pdf",bbox_inches="tight",pad_inches=0)

                # ax.set_xlim(40)
                # ax.set_ylim(40)
                # for r in robots:
                #     if search_state:
                #         intruder.body.set_visible(False)
                #         for rb in robots:
                #             rb.body.set_visible(False)
                #             if rb.plott!= None:
                #                 rb.plott[0].set_visible(False)
                # plt.show()
                # plt.pause(0.0000001)    
            data_runs.append([num_robots,t])
        data_his.append(data_runs)
    bag = data_his
    fileObject = open(path_+'/results_7s_5g_10lg/0_MRIS_s_bench', 'wb')
    pkl.dump(bag,fileObject)
    fileObject.close()

dynamic_intruder = 1
if dynamic_intruder:
    path_ = os.path.realpath(os.path.dirname(__file__))
    # fig,ax = plt.subplots()
    # plt.box(False)
    # ax.set_aspect('equal')
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    ax = None
    itera = np.random.randint(50)
    print('itera :',itera)
    tile_sixe = 10
    grid_height = 5
    grid_width = 5
    
    # Random simple connected rectilinear polygon
    boundary,edges,start_edge,end_edge = Inflate_Cut_algorithm(20,g_size=tile_sixe)#,ax=ax)
    boundary1 = get_area_vertices(np.concatenate((boundary[2:],boundary[:3]),axis=0))
    # ax.clear()
    # ax.plot(boundary1[:,0],boundary1[:,1],color='indigo')
    
    grid_graph = Grid_graph(boundary1,grid_size=grid_height,ax=ax)
    
    # plt.ioff()
    # plt.show()

    data_his = []
    # plt.ion()
    grid_mesh = None
    for num_robots in range(1,153,3):#range(1,39):#range(1,153,3):
        avg_time = []
        data_runs = []
        print(num_robots)
        for runs in range(100):
            for g in grid_graph:
                grid_graph[g].robo_path_pre = None
                grid_graph[g].heuristics = 0
                grid_graph[g].distance = 0
                grid_graph[g].total_cost = 0
            t = 0   #   Initialize time
            search_state =  0 #   0 = not found, 1 = found


            robots = searchers(num_robots,grid_width,grid_height,grid_graph,ax)
            intruder = intruder_(grid_graph,grid_width,grid_height,ax)

            while not search_state: #   While Intruder is not found
                for r in robots:
                    if intruder.present_loc == r.present_loc:
                        # intruder.body.set_visible(False)
                        # if intruder.plott!= None:
                        #     intruder.plott[0].set_visible(False)
                        # for rb in robots:
                        #     rb.body.set_visible(False)
                        #     if rb.plott!= None:
                        #         rb.plott[0].set_visible(False)
                        search_state = 1
                        avg_time.append(t)
                present_locs = [r.present_loc for r in robots]
                update_prob_costs(grid_graph,present_locs)
                # if grid_mesh is None:
                #     grid_mesh = plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=0)
                # else:
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=1)
                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].next_loc = intruder.present_loc
                        for g in grid_graph:
                            grid_graph[g].robo_path_pre = None
                            
                        if robots[i].next_loc != robots[i].present_loc:
                            path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)
                        else:
                            path = [robots[i].next_loc]

                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_ahead = path
                        robots[i].path_progressor = 0
                    

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    # robots[i].update(grid_graph)
                    robots[i].path_progressor += 1
                    if len(robots[i].trajectory)<5:
                        robots[i].trajectory.append([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                        # if len(robots[i].trajectory)==4:
                        #     robots[i].plott = ax.plot(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                    else:
                        robots[i].update_trajectory([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                        # robots[i].plott[0].set_data(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1])

                # Intruder update
                intruder.next_loc = intruder.loc(grid_graph)
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc
                # intruder.update(grid_graph)
                if len(intruder.trajectory)<10:
                    intruder.trajectory.append([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                    # if len(intruder.trajectory)==9:
                    #     intruder.plott = ax.plot(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1],color='red',linewidth=3,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))
                else:
                    intruder.update_trajectory([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                    # intruder.plott[0].set_data(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1])
                t += 1
                # if runs==1 and t==60:
                #     path = os.getcwd()
                #     plt.savefig(path+'/results/1RIS_d_random.pdf',format = "pdf",bbox_inches="tight",pad_inches=0)

                # ax.set_xlim(40)
                # ax.set_ylim(40)
                # for r in robots:
                #     if search_state:
                #         intruder.body.set_visible(False)
                #         if intruder.plott!= None:
                #             intruder.plott[0].set_visible(False)
                #         for rb in robots:
                #             rb.body.set_visible(False)
                #             if rb.plott!= None:
                #                 rb.plott[0].set_visible(False)
                # plt.show()
                # plt.pause(0.000001)     
            data_runs.append([num_robots,t])
        data_his.append(data_runs)
    bag = data_his
    fileObject = open(path_+'/results_7s_5g_10lg/0_MRIS_d_bench', 'wb')
    pkl.dump(bag,fileObject)
    fileObject.close()
