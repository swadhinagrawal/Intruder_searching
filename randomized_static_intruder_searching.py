#   Authors: Swadhin Agrawal

import matplotlib.pyplot as plt
import numpy as np
import copy as cp
import pickle as pkl
import os
np.random.seed(748657)

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

def Inflate_Cut_algorithm(num_vertices,ax=None):
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

    def Inflate(p,c,ax=None):
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
            p_trial = Inflate(p_trial,random_c)

            p_trial, cut_success = Cut(p_trial,random_c)
        P = p_trial
        ax1.clear()
        polygon_boundary(P,ax1)
        r -= 1
    boundary,edges,start_edge,end_edge = polygon_boundary(P,ax)
    return boundary,edges,start_edge,end_edge , P

# Decomposition

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

def get_rectangles(p,vertex_set,ax):
    grids = cp.copy(p)
    rectangles = []
    while len(grids)>0:
        g = list(grids.items())[0][0]
        marker = [0,0,0,0]  # l,r,t,b
        if grids[g].l_nei is not None:
            marker[0] = 1
        if grids[g].r_nei is not None:
            marker[1] = 1
        if grids[g].t_nei is not None:
            marker[2] = 1
        if grids[g].b_nei is not None:
            marker[3] = 1
        while grids[g].l_nei != None:
            if grids[g].l_nei != None:
                g = grids[g].l_nei
            else:
                break
        while grids[g].b_nei != None:
            if grids[g].b_nei != None:
                g = grids[g].b_nei
            else:
                break

        # x = grids[g].get_x()
        # y = grids[g].get_y()
        # ax.plot(x,y,color='red')
        plt.show()
        def check_presence_inside(vm,vertex_set,C):
            a = False
            for p2 in vertex_set:
                a = a or (p2[0]<vm[0] and p2[1]<vm[1] and p2[0]>C[0] and p2[1]>C[1])
            return a

        rec_vertices = [grids[g].bl,grids[g].tl,grids[g].br]
        min_x = min(np.array(rec_vertices)[:,0])
        min_y = min(np.array(rec_vertices)[:,1])
        max_x = max(np.array(rec_vertices)[:,0])
        max_y = max(np.array(rec_vertices)[:,1])
        g_h = g
        g_v = g
        g_v1 = g_v
        g_h1 = g_h
        while grids[g_h].r_nei!= None or grids[g_v].t_nei!= None:
            changed = 0
            if np.sum(marker)>0:
                
                if grids[g_h].r_nei!= None:
                    
                    g_h = grids[g_h].r_nei
                    
                    rec_vertices.append(grids[g_h].br)
                    max_x = max(np.array(rec_vertices)[:,0])
                    max_y = max(np.array(rec_vertices)[:,1])
                    if not check_presence_inside(np.array([max_x,max_y]),vertex_set,np.array([min_x,min_y])):
                        g_h1 = g_h
                        changed = 1

                    # x = grids[g_h].get_x()
                    # y = grids[g_h].get_y()
                    # ax.plot(x,y,color='pink')
                    # plt.show()
                    
                if grids[g_v].t_nei!= None:
                    
                    g_v = grids[g_v].t_nei
                    rec_vertices.append(grids[g_v].tl)
                    max_x = max(np.array(rec_vertices)[:,0])
                    max_y = max(np.array(rec_vertices)[:,1])
                    if not check_presence_inside(np.array([max_x,max_y]),vertex_set,np.array([min_x,min_y])):
                        g_v1 = g_v
                        changed = 1
                    # x = grids[g_v].get_x()
                    # y = grids[g_v].get_y()
                    # ax.plot(x,y,color='pink')
                    # plt.show()
                    
            else:
                rec_vertices = grids[g].get_corners()
                rec_vertices = np.concatenate((rec_vertices,np.array([rec_vertices[0]])),axis=0)
            if not changed:
                break

        max_x = grids[g_h1].br[0]

        max_y = grids[g_v1].tl[1]

        rectangles.append(Grid(tl=np.array([min_x,max_y]),bl=np.array([min_x,min_y]),tr=np.array([max_x,max_y]),br=np.array([max_x,min_y])))
        if ax != None:
            ax.plot(rectangles[-1].get_x(),rectangles[-1].get_y())
        delete_these = []
        delete_border_grids = []
        for g in grids:
            x = grids[g].get_x()
            y = grids[g].get_y()
            marker = [0,0,0,0]  # l,r,t,b
            if grids[g].l_nei is not None:
                marker[0] = 1
            if grids[g].r_nei is not None:
                marker[1] = 1
            if grids[g].t_nei is not None:
                marker[2] = 1
            if grids[g].b_nei is not None:
                marker[3] = 1
            if min(x) > min(rectangles[-1].get_x()) and max(x) < max(rectangles[-1].get_x()) and min(y) > min(rectangles[-1].get_y()) and max(y) < max(rectangles[-1].get_y()) and np.sum(marker)==4:
                delete_these.append(g)
            elif min(x) >= min(rectangles[-1].get_x()) and max(x) <= max(rectangles[-1].get_x()) and min(y) >= min(rectangles[-1].get_y()) and max(y) <= max(rectangles[-1].get_y()):
                delete_border_grids.append(g)
                p[g].rectangle_num = len(rectangles)-1
        
        for gg in delete_these:
            if grids[gg].t_nei != None:
                grids[grids[gg].t_nei].b_nei = None
            if grids[gg].b_nei != None:
                grids[grids[gg].b_nei].t_nei = None
            if grids[gg].r_nei != None:
                grids[grids[gg].r_nei].l_nei = None
            if grids[gg].l_nei != None:
                grids[grids[gg].l_nei].r_nei = None
        for gg in delete_these:
            del grids[gg]
        plt.show()
        
        def get_neighbour(cel,s):
            if s == 'r':
                return cel.r_nei
            elif s == 'l':
                return cel.l_nei
            elif s == 't':
                return cel.t_nei
            elif s == 'b':
                return cel.b_nei
        
        sides = ['l','r','t','b']
        if ax!=None:
            ax.clear()
            ax.plot(vertex_set[:,0],vertex_set[:,1],color='darkblue')
        for gg in delete_border_grids:
            for s in sides:
                nei = get_neighbour(p[gg],s)
                if nei not in delete_border_grids and nei != None:
                    if grids[nei].t_nei in delete_border_grids:
                        p[nei].passage = 1
                        p[nei].passage_t = 1
                        if ax!=None:
                            ax.plot(p[nei].get_x()[2:-1],p[nei].get_y()[2:-1],color='black')
                    if grids[nei].b_nei in delete_border_grids:
                        p[nei].passage = 1
                        p[nei].passage_b = 1
                        if ax!=None:
                            ax.plot(p[nei].get_x()[:2],p[nei].get_y()[:2],color='black')
                    if grids[nei].l_nei in delete_border_grids:
                        p[nei].passage = 1
                        p[nei].passage_l = 1
                        if ax!=None:
                            ax.plot([p[nei].get_x()[ind] for ind in range(-2,1)],[p[nei].get_y()[ind] for ind in range(-2,1)],color='black')
                    if grids[nei].r_nei in delete_border_grids:
                        p[nei].passage = 1
                        p[nei].passage_r = 1
                        if ax!=None:
                            ax.plot(p[nei].get_x()[1:3],p[nei].get_y()[1:3],color='black')
        
        for gg in delete_border_grids:
            if grids[gg].t_nei != None:
                grids[grids[gg].t_nei].b_nei = None
            if grids[gg].b_nei != None:
                grids[grids[gg].b_nei].t_nei = None
            if grids[gg].r_nei != None:
                grids[grids[gg].r_nei].l_nei = None
            if grids[gg].l_nei != None:
                grids[grids[gg].l_nei].r_nei = None
                
        for gg in delete_border_grids:
            del grids[gg]

        # plt.show()
    
    # for g in p:
    #     # ax.plot(p[g].get_x(),p[g].get_y(),color='orange')
    #     if p[g].passage_l:
    #         ax.plot([p[g].get_x()[ind] for ind in range(-2,1)],[p[g].get_y()[ind] for ind in range(-2,1)],color='white')
    #     if p[g].passage_r:
    #         ax.plot(p[g].get_x()[1:3],p[g].get_y()[1:3],color='white')
    #     if p[g].passage_t:
    #         ax.plot(p[g].get_x()[2:-1],p[g].get_y()[2:-1],color='white')
    #     if p[g].passage_b:
    #         ax.plot(p[g].get_x()[:2],p[g].get_y()[:2],color='white')
    # plt.show()
    return rectangles, p

#   Spacefilling curves
def make_grid(x1,y1,x2,y2,inner_grid_height,inner_grid_width,ax,grid_graph,grids,rec_num):
    def check_presence_inside(vm,vertex_set,C):
        a = False
        for p2 in vertex_set:
            a = a or (p2[0]<vm[0] and p2[1]<vm[1] and p2[0]>C[0] and p2[1]>C[1])
        return a
    centroids = []
    for i in range(x1,x2,inner_grid_width):
        for j in range(y1,y2,inner_grid_height):
            grid_graph.update({str([i,j]):Grid(tl=np.array([i,j+inner_grid_height]),tr=np.array([i+inner_grid_width,j+inner_grid_height]),bl=np.array([i,j]),br=np.array([i+inner_grid_width,j]))})
            grid_graph[str([i,j])].height = inner_grid_height
            grid_graph[str([i,j])].width = inner_grid_width
            centroids.append(grid_graph[str([i,j])].centroid)
            if i-inner_grid_width>=x1:
                grid_graph[str([i,j])].l_nei = str([i-inner_grid_width,j])
            if i+inner_grid_width<x2:
                grid_graph[str([i,j])].r_nei = str([i+inner_grid_width,j])
            if j-inner_grid_height>=y1:
                grid_graph[str([i,j])].b_nei = str([i,j-inner_grid_height])
            if j+inner_grid_height<y2:
                grid_graph[str([i,j])].t_nei = str([i,j+inner_grid_height])
            if i == 99 and j == 90:
                print(' ')
            if i == x1 or i == x2-inner_grid_width or j == y1 or j == y2-inner_grid_height:
                for g in grids:
                    if check_presence_inside([x2,y2],[grids[g].centroid],[x1,y1]) and (grids[g].passage == 1) and check_presence_inside(grids[g].tr,[grid_graph[str([i,j])].centroid],grids[g].bl):
                        if grids[g].passage_l and grid_graph[str([i,j])].centroid[0] == grids[g].bl[0]+inner_grid_width/2:
                            grid_graph[str([i,j])].passage = 1
                            grid_graph[str([i,j])].l_nei = str(list((grid_graph[str([i,j])].bl-np.array([inner_grid_width,0])).astype(int)))
                            grid_graph[grid_graph[str([i,j])].l_nei].r_nei = str([i,j])
                        if grids[g].passage_r and grid_graph[str([i,j])].centroid[0] == grids[g].br[0]-inner_grid_width/2:
                            grid_graph[str([i,j])].passage = 1
                            grid_graph[str([i,j])].r_nei = str(list((grid_graph[str([i,j])].br).astype(int)))
                            grid_graph[grid_graph[str([i,j])].r_nei].l_nei = str([i,j])
                        if grids[g].passage_t and grid_graph[str([i,j])].centroid[1] == grids[g].tl[1]-inner_grid_height/2:
                            grid_graph[str([i,j])].passage = 1
                            grid_graph[str([i,j])].t_nei = str(list((grid_graph[str([i,j])].tl).astype(int)))
                            grid_graph[grid_graph[str([i,j])].t_nei].b_nei = str([i,j])
                        if grids[g].passage_b and grid_graph[str([i,j])].centroid[1] == grids[g].bl[1]+inner_grid_height/2:
                            grid_graph[str([i,j])].passage = 1
                            grid_graph[str([i,j])].b_nei = str(list((grid_graph[str([i,j])].bl-np.array([0,inner_grid_height])).astype(int)))
                            grid_graph[grid_graph[str([i,j])].b_nei].t_nei = str([i,j])

    return grid_graph,centroids

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
                temp_grids[child].total_cost = temp_grids[current_node].total_cost + temp_grids[child].heuristics# + 0.1*temp_grids[child].heuristics_2

            # Add the child to the open list
            if child not in open_nodes+closed_nodes and child != None:
                open_nodes.append(child)

    return path

# def plot_mesh(fig,ax,grid,title,without_bar=0):
#     values_x = []
#     values_y = []
#     for g in grid:
#         values_x = values_x + list(grid[g].get_x())
#         values_y = values_y + list(grid[g].get_y())

#     x = np.arange(min(values_x),max(values_x),grid_width)
#     y = np.arange(min(values_y),max(values_y),grid_height)
#     z = []
#     for i in x:
#         verti = []
#         for j in y:
#             if str([i, j]) in grid:
#                 verti.append(grid[str([i,j])].heuristics)
#             else:
#                 verti.append(0)
#         z.append(verti)
#     z = np.array(z).T
#     cs = ax.pcolormesh(x + grid_width/2,y + grid_height/2,z,shading='auto',cmap='inferno')
#     cs.set_clim(0,1)
#     if not without_bar:
#         cbar = fig.colorbar(cs,orientation='vertical')
#         cbar.set_ticks(np.arange(0,1.1,0.1))

#     ax.set_aspect('equal', 'box')
#     # ax.set_xticks([])
#     # ax.set_yticks([])
#     ax.set_title(title)

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
            if str([i, j]) in grid:
                verti.append(grid[str([i,j])].heuristics)
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
    for g in grids: #   Update the Costs and Probabilities of the grids
        for l in present_locs:
            if l==g:
                # grids[g].probability -= 0#1/(arena_h*arena_w)
                # if grids[g].probability<0:
                #     grids[g].probability=0
                # elif grids[g].probability>1:
                #     grids[g].probability=1
                grids[g].heuristics += 50/len(grids)
            # else:
            #     grids[g].probability += 0#1/(arena_h*arena_w*len(present_locs))

# Single robot
_1_static_intruder_random_arena = 0
if _1_static_intruder_random_arena:
    path_ = os.getcwd()
    # fig,ax = plt.subplots()
    # plt.box(False)
    # ax.set_aspect('equal')
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(np.random.randint(30))#,ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        ax = None
        grid_graph,centroids = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)


    time = []
    agents = []
    # plt.ion()
    grid_mesh = None
    once = 0
    for num_robots in range(1,2):#10,1):
        avg_time = []
        for runs in range(10):
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
                        # for rb in robots:
                        #     rb.body.set_visible(False)
                        #     if rb.plott!= None:
                        #         rb.plott[0].set_visible(False)
                        search_state = 1
                        avg_time.append(t)
                present_locs = [r.present_loc for r in robots]
                update_prob_costs(grid_graph,present_locs)
                # if not once:
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=0)
                #     once = 1
                # else:    
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=1)

                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].next_loc = robots[i].loc(grid_graph)
                        while robots[i].next_loc == robots[i].present_loc:
                            robots[i].next_loc = robots[i].loc(grid_graph)
                        for g in grid_graph:
                            grid_graph[g].robo_path_pre = None
                        path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)

                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_ahead = path
                        robots[i].path_progressor = 0
                    # robots[i].update(grid_graph)

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    robots[i].path_progressor += 1

                    # if len(robots[i].trajectory)<5:
                    #     robots[i].trajectory.append([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     if len(robots[i].trajectory)==4:
                    #         robots[i].plott = ax.plot(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                    # else:
                    #     robots[i].update_trajectory([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     robots[i].plott[0].set_data(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1])

                # Intruder update
                intruder.next_loc = intruder.present_loc
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc
                t += 1
                
                # ax.set_xlim(40,190)
                # ax.set_ylim(40,190)
                # plt.show()
                # plt.pause(0.000001)     
            
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

    fileObject = open(path_+'/results/1RIS_random', 'wb')
    pkl.dump(np.array(avg_time),fileObject)

_1_dynamic_intruder_random_arena = 0
if _1_dynamic_intruder_random_arena:
    path_ = os.getcwd()
    # fig,ax = plt.subplots()
    # plt.box(False)
    # ax.set_aspect('equal')
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(np.random.randint(30))#,ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        ax = None
        grid_graph,centroids = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)


    time = []
    agents = []
    # plt.ion()
    once = 0
    grid_mesh = None
    for num_robots in range(1,2):
        avg_time = []
        for runs in range(10):
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
                # if not once:
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=0)
                #     once = 1
                # else:
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=1)
                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].next_loc = robots[i].loc(grid_graph)
                        while robots[i].next_loc == robots[i].present_loc:
                            robots[i].next_loc = robots[i].loc(grid_graph)
                        for g in grid_graph:
                            grid_graph[g].robo_path_pre = None
                        path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)

                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_ahead = path
                        robots[i].path_progressor = 0
                    # robots[i].update(grid_graph)

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    robots[i].path_progressor += 1
                    # if len(robots[i].trajectory)<5:
                    #     robots[i].trajectory.append([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     if len(robots[i].trajectory)==4:
                    #         robots[i].plott = ax.plot(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                    # else:
                    #     robots[i].update_trajectory([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     robots[i].plott[0].set_data(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1])

                # Intruder update
                intruder.next_loc = intruder.loc(grid_graph)
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc
                # intruder.update(grid_graph)
                # intruder.body.set_offsets([[grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]]])
                # if len(intruder.trajectory)<10:
                #     intruder.trajectory.append([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                #     if len(intruder.trajectory)==9:
                #         intruder.plott = ax.plot(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1],color='red',linewidth=3,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))
                # else:
                #     intruder.update_trajectory([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                #     intruder.plott[0].set_data(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1])

                t += 1

                # ax.set_xlim(40,190)
                # ax.set_ylim(40,190)
                # plt.show()
                # plt.pause(0.000001)     
            
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

    fileObject = open(path_+'/results/1RIS_random_d', 'wb')
    pkl.dump(np.array(avg_time),fileObject)
    fileObject.close()



# Multi robot
static_intruder_random_arena = 0
if static_intruder_random_arena:
    path_ = os.getcwd()
    # fig,ax = plt.subplots()
    # plt.box(False)
    # ax.set_aspect('equal')
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(np.random.randint(30))#,ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        ax = None
        grid_graph,centroids = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)


    time = []
    agents = []
    # plt.ion()
    grid_mesh = None
    once = 0
    for num_robots in range(1,10,1):
        avg_time = []
        for runs in range(10):
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
                        # for rb in robots:
                        #     rb.body.set_visible(False)
                        #     if rb.plott!= None:
                        #         rb.plott[0].set_visible(False)
                        search_state = 1
                        avg_time.append(t)
                present_locs = [r.present_loc for r in robots]
                update_prob_costs(grid_graph,present_locs)
                # if not once:
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=0)
                #     once = 1
                # else:    
                #     plot_mesh(fig,ax,grid_graph,'cost map',grid_mesh,without_bar=1)

                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].next_loc = robots[i].loc(grid_graph)
                        while robots[i].next_loc == robots[i].present_loc:
                            robots[i].next_loc = robots[i].loc(grid_graph)
                        for g in grid_graph:
                            grid_graph[g].robo_path_pre = None
                        path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)

                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_ahead = path
                        robots[i].path_progressor = 0
                    # robots[i].update(grid_graph)

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    robots[i].path_progressor += 1

                    # if len(robots[i].trajectory)<5:
                    #     robots[i].trajectory.append([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     if len(robots[i].trajectory)==4:
                    #         robots[i].plott = ax.plot(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                    # else:
                    #     robots[i].update_trajectory([grid_graph[robots[i].present_loc].centroid[0],grid_graph[robots[i].present_loc].centroid[1]])
                    #     robots[i].plott[0].set_data(np.array(robots[i].trajectory)[:,0],np.array(robots[i].trajectory)[:,1])

                # Intruder update
                intruder.next_loc = intruder.present_loc
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc
                t += 1
                
                # ax.set_xlim(40,190)
                # ax.set_ylim(40,190)
                # plt.show()
                # plt.pause(0.000001)     
            
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

    fileObject = open(path_+'/results/MRIS_random_s', 'wb')
    pkl.dump(np.array(time),fileObject)
    pkl.dump(np.array(agents),fileObject)
    fileObject.close()

dynamic_intruder_random_arena = 0
if dynamic_intruder_random_arena:
    path_ = os.getcwd()
    # fig,ax = plt.subplots()
    # ax.set_aspect('equal')
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(np.random.randint(30))#,ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        ax = None
        grid_graph,centroids = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)


    time = []
    agents = []
    plt.ion()

    for num_robots in range(1,10,1):
        avg_time = []
        for runs in range(10):
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

                present_locs = [r.present_loc for r in robots]
                update_prob_costs(grid_graph,present_locs)
                # if num_robots == 5 and t == 20:
                #     plot_mesh(fig,ax,grid_graph,'cost map',without_bar=0)

                #   Robots update
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if robots[i].present_loc == robots[i].path_ahead[-1]:
                        robots[i].next_loc = robots[i].loc(grid_graph)
                        while robots[i].next_loc == robots[i].present_loc:
                            robots[i].next_loc = robots[i].loc(grid_graph)
                        for g in grid_graph:
                            grid_graph[g].robo_path_pre = None
                        path = Astar(grid_graph,robots[i].present_loc,robots[i].next_loc)#,fig,ax)

                        robots[i].path_history = robots[i].path_history + robots[i].path_ahead
                        robots[i].path_ahead = path
                        robots[i].path_progressor = 0
                    # robots[i].update(grid_graph)

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_ahead[robots[i].path_progressor]
                    robots[i].path_progressor += 1
                # Intruder update
                intruder.next_loc = intruder.loc(grid_graph)
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc
                # intruder.update(grid_graph)
                t += 1
                for r in robots:
                    if intruder.present_loc == r.present_loc:
                        # intruder.body.set_visible(False)
                        # for rb in robots:
                        #     rb.body.set_visible(False)
                        search_state = 1
                        avg_time.append(t)
                # plt.show()
                # plt.pause(0.000001)     
            
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

    fileObject = open(path_+'/results/MRIS_random_d', 'wb')
    pkl.dump(np.array(time),fileObject)
    pkl.dump(np.array(agents),fileObject)
    fileObject.close()


# Results plotter
single_searcher_plot_results = 0
if single_searcher_plot_results:
    path = os.getcwd()
    file = open(path+'/results/MRIS_random_s', 'rb')
    obj = pkl.load(file)

    plt.ioff()
    fig,ax = plt.subplots()
    
    ax.plot(obj,range(1,10),linewidth=4)
    # ax.plot(range(10),obj,linewidth=4)
    # ax.plot(range(9),[np.sum(obj)/len(obj) for i in range(10)],linewidth=4)
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    plt.title('SRIS: Random; static intruder',fontdict=t_font)
    plt.xlabel('Runs',fontdict=ax_font)
    plt.ylabel('Search Time',fontdict=ax_font)
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    plt.savefig(path+'/results/MRIS_random_s_plot.png')
    plt.show()

multi_searcher_plot_results = 1
if multi_searcher_plot_results:
    path = os.getcwd()
    file = open(path+'/results/MRIS_random_s', 'rb')
    obj = pkl.load(file)

    plt.ioff()
    fig,ax = plt.subplots()
    
    ax.plot(obj,range(1,10),linewidth=4,label='Static Intruder')
    # ax.plot(range(10),obj,linewidth=4)
    # ax.plot(range(9),[np.sum(obj)/len(obj) for i in range(10)],linewidth=4)
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    plt.title('SRIS: Random; static intruder',fontdict=t_font)
    plt.xlabel('Runs',fontdict=ax_font)
    plt.ylabel('Search Time',fontdict=ax_font)
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    file.close()
    # plt.savefig(path+'/results/MRIS_random_s_plot.png')
    # plt.show()
    file = open(path+'/results/MRIS_random_d', 'rb')
    obj = pkl.load(file)

    plt.ioff()
    # fig,ax = plt.subplots()
    
    ax.plot(obj,range(1,10),linewidth=4,label='Dynamic Intruder')
    # ax.plot(range(10),obj,linewidth=4)
    # ax.plot(range(9),[np.sum(obj)/len(obj) for i in range(10)],linewidth=4)
    t_font = {'weight': 'bold',
        'size': 15}
    ax_font = {'weight': 'bold',
        'size': 20}
    plt.title('SRIS: Random; static intruder',fontdict=t_font)
    plt.xlabel('Runs',fontdict=ax_font)
    plt.ylabel('Search Time',fontdict=ax_font)
    plt.xticks(fontsize=15,fontweight='bold')
    plt.yticks(fontsize=15,fontweight='bold')
    plt.tight_layout()
    plt.legend()
    plt.savefig(path+'/results/MRIS_random_sd_plot.png')
    plt.show()

# static_intruder = 0
# if static_intruder:
#     arena_w = 10
#     arena_h = 10

#     time = []
#     agents = []
#     plt.ion()
#     fig,ax = plt.subplots()
#     fig1,ax1 = plt.subplots()
#     faulthandler.enable()
#     for num_robots in range(1,25,1):
#         avg_time = []
#         for runs in range(10):
#             t = 0   #   Initialize time
#             search_state =  0 #   0 = not found, 1 = found
#             ax.clear()
#             ax1.clear()
#             grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))
#             probabilities = probabilities_(grids)
#             costs = cost_map_(grids)
#             if num_robots ==1 and runs ==0:
#                 without_bar = 0
#             else:
#                 without_bar = 1
            
#             plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar)
#             plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar)    

#             robots = searchers(num_robots,arena_w,arena_h,cp.copy(grids),ax)
#             intruder = intruder_(arena_w,arena_h,cp.copy(grids),ax)
#             intruder_plot(intruder,ax1)
#             searchers_plot(robots,ax1)
            
#             while not search_state: #   While Intruder is not found
#                 ax.clear()
#                 ax1.clear()
#                 present_locs = visited_loc(robots)
#                 update_prob_costs(grids,present_locs,arena_h,arena_w)
#                 probabilities = probabilities_(grids)
#                 costs = cost_map_(grids)
#                 plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar=1)
#                 plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar=1)

#                 intruder_plot(intruder,ax1)
#                 searchers_plot(robots,ax1)

#                 plot_past_path(ax,robots)
#                 plt.show()
#                 plt.pause(0.00001)

#                 # Intruder update
#                 intruder.next_loc = intruder.present_loc

#                 intruder.past_loc = intruder.present_loc
#                 intruder.present_loc = intruder.next_loc

#                 ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],edgecolors='white',s=25, facecolors='none')
#                 ax.scatter([intruder.next_loc[1]+0.5],[intruder.next_loc[0]+0.5],marker='x',color = 'green')
#                 plt.show()
#                 plt.pause(0.00001)
                

#                 #   Robots update

#                 # paths = []
#                 for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
#                     if np.linalg.norm(robots[i].present_loc-robots[i].path_trace[-1])==0:
#                         robots[i].next_loc = robots[i].loc(grids)

#                         for g in grids:
#                             for g_ij in g:
#                                 g_ij.parent = None
#                         path = Astar(grids,robots[i].present_loc,robots[i].next_loc)
#                         for p in range(len(path)-1):
#                             ax.scatter([path[p][1]+0.5],[path[p][0]+0.5],edgecolors='cyan',s=25, facecolors='none')

#                             # ax.plot([path[p][1]+0.5,path[p+1][1]+0.5],[path[p][0]+0.5,path[p+1][0]+0.5],color= 'black',alpha=0.2)
#                         ax.scatter([path[0][1]+0.5],[path[0][0]+0.5],edgecolors='white',s=25, facecolors='none')
#                         ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],marker='x',color = 'green')
#                         plt.show()
#                         plt.pause(0.00001)
                        
#                         robots[i].path_trace_past = path#robots[i].path_trace
#                         robots[i].path_trace = path
#                         robots[i].path_progressor = 0

#                     # paths+= list(robots[i].present_loc)

#                     robots[i].past_loc = robots[i].present_loc
                    
#                     robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment
#                     robots[i].path_progressor += 1

#                 t += 1
#                 for r in robots:
#                     if np.linalg.norm(intruder.present_loc-r.present_loc)==0:
#                         search_state = 1
#                         break  
#                 plt.show()
#                 plt.pause(0.000001)     
#             avg_time.append(t)
#         time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
#         agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

#     fileObject = open('data_static_intruder', 'wb')
#     pkl.dump(np.array(time),fileObject)
#     pkl.dump(np.array(agents),fileObject)
#     plt.ioff()
#     fig,ax = plt.subplots()
#     ax.plot(time,agents)
#     ax.set_xlabel("Time taken")
#     ax.set_ylabel("Number of search agents")
#     plt.savefig('static_intruder.eps')
#     plt.show()

# dynamic_intruder = 0
# if dynamic_intruder:
#     arena_w = 10
#     arena_h = 10

#     time = []
#     agents = []
#     plt.ion()
#     fig,ax = plt.subplots()
#     fig1,ax1 = plt.subplots()
#     faulthandler.enable()
#     for num_robots in range(1,10,1):
#         avg_time = []
#         for runs in range(10):
#             t = 0   #   Initialize time
#             search_state =  0 #   0 = not found, 1 = found
#             ax.clear()
#             ax1.clear()
#             grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))
#             probabilities = probabilities_(grids)
#             costs = cost_map_(grids)
#             if num_robots ==1 and runs ==0:
#                 without_bar = 0
#             else:
#                 without_bar = 1
            
#             plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar)
#             plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar)    

#             robots = searchers(num_robots,arena_w,arena_h,cp.copy(grids),ax)
#             intruder = intruder_(arena_w,arena_h,cp.copy(grids),ax)
#             intruder_plot(intruder,ax1)
#             searchers_plot(robots,ax1)
            
#             while not search_state: #   While Intruder is not found
#                 ax.clear()
#                 ax1.clear()
#                 present_locs = visited_loc(robots)
#                 update_prob_costs(grids,present_locs,arena_h,arena_w)
#                 probabilities = probabilities_(grids)
#                 costs = cost_map_(grids)
#                 plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar=1)
#                 plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar=1)

#                 intruder_plot(intruder,ax1)
#                 searchers_plot(robots,ax1)

#                 plot_past_path(ax,robots)
#                 plt.show()
#                 plt.pause(0.00001)

#                 # Intruder update
#                 intruder.next_loc = intruder.loc(grids)

#                 intruder.past_loc = intruder.present_loc
#                 intruder.present_loc = intruder.next_loc

#                 ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],edgecolors='white',s=25, facecolors='none')
#                 ax.scatter([intruder.next_loc[1]+0.5],[intruder.next_loc[0]+0.5],marker='x',color = 'green')
#                 plt.show()
#                 plt.pause(0.00001)
                

#                 #   Robots update

#                 paths = []
#                 for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
#                     if np.linalg.norm(robots[i].present_loc-robots[i].path_trace[-1])==0:
#                         robots[i].next_loc = robots[i].loc(grids)

#                         for g in grids:
#                             for g_ij in g:
#                                 g_ij.parent = None
#                         path = Astar(grids,robots[i].present_loc,robots[i].next_loc)
#                         for p in range(len(path)-1):
#                             ax.scatter([path[p][1]+0.5],[path[p][0]+0.5],edgecolors='cyan',s=25, facecolors='none')

#                             # ax.plot([path[p][1]+0.5,path[p+1][1]+0.5],[path[p][0]+0.5,path[p+1][0]+0.5],color= 'black',alpha=0.2)
#                         ax.scatter([path[0][1]+0.5],[path[0][0]+0.5],edgecolors='white',s=25, facecolors='none')
#                         ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],marker='x',color = 'green')
#                         plt.show()
#                         plt.pause(0.00001)
                        
#                         robots[i].path_trace_past = path#robots[i].path_trace
#                         robots[i].path_trace = path
#                         robots[i].path_progressor = 0

#                     paths+= list(robots[i].present_loc)

#                     robots[i].past_loc = robots[i].present_loc

#                     robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment
#                     robots[i].path_progressor += 1

#                 t += 1
#                 for r in robots:
#                     if np.linalg.norm(intruder.present_loc-r.present_loc)==0:
#                         search_state = 1
#                         break  
#                 plt.show()
#                 plt.pause(0.000001)     
#             avg_time.append(t)
#         time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
#         agents.append(num_robots)   #   Store number of robots utilized for search of static intruder
    
#     fileObject = open('data_dynamic_intruder', 'wb')
#     pkl.dump(np.array(time),fileObject)
#     pkl.dump(np.array(agents),fileObject)
#     plt.ioff()
#     fig,ax = plt.subplots()
#     ax.plot(time,agents)
#     ax.set_xlabel("Time taken")
#     ax.set_ylabel("Number of search agents")
#     plt.savefig('dynamic_intruder.eps')
#     plt.show()
