#   Authors: Aayush Gohil, Swadhin Agrawal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy as cp
import pickle as pkl
import os
from random import randint,random,sample
from shapely.geometry import LineString

sed = 19778167 # 5 spike 19778167: 38lg
# sed = 27819661 # 7 spike 27819661 38lg
np.random.seed(sed)

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
        
        self.path_pre = None
        self.path_post = None
        self.centroid = (self.tl+self.bl+self.tr+self.br)/4.0

        self.robot_home = False
        self.passage = 0
        self.passage_l = 0
        self.passage_r = 0
        self.passage_t = 0
        self.passage_b = 0

        self.rectangle_num = None
        self.border_grid = 0
        self.path_connector = None
    
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
    return boundary,edges,start_edge,end_edge , P

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
            tl = np.array([i,j+grid_size])
            tr = np.array([i+grid_size,j+grid_size])
            bl = np.array([i,j])
            br = np.array([i+grid_size,j])
            g = Grid(tl=tl,bl=bl,tr=tr,br=br)
            if is_grid_in_poly(polygon,g,min_x,min_y,max_x,max_y):
                if str(bl - np.array([grid_size,0])) in grids:
                    g.l_nei = str(bl - np.array([grid_size,0]))
                    grids[str(bl - np.array([grid_size,0]))].r_nei = str(bl)
                    ax.plot([g.centroid[0],grids[str(bl - np.array([grid_size,0]))].centroid[0]],[g.centroid[1],grids[str(bl - np.array([grid_size,0]))].centroid[1]],color='maroon',linewidth=0.5)
                if str(br) in grids:
                    g.r_nei = str(br)
                    grids[str(br)].l_nei = str(bl)
                    ax.plot([g.centroid[0],grids[str(br)].centroid[0]],[g.centroid[1],grids[str(br)].centroid[1]],color='maroon',linewidth=0.5)
                if str(tl) in grids:
                    g.t_nei = str(tl)
                    grids[str(tl)].b_nei = str(bl)
                    ax.plot([g.centroid[0],grids[str(tl)].centroid[0]],[g.centroid[1],grids[str(tl)].centroid[1]],color='maroon',linewidth=0.5)
                if str(bl - np.array([0,grid_size])) in grids:
                    g.b_nei = str(bl - np.array([0,grid_size]))
                    grids[str(bl - np.array([0,grid_size]))].t_nei = str(bl)
                    ax.plot([g.centroid[0],grids[str(bl - np.array([0,grid_size]))].centroid[0]],[g.centroid[1],grids[str(bl - np.array([0,grid_size]))].centroid[1]],color='maroon',linewidth=0.5)
                grids[str(bl)] = g
                x = grids[str(bl)].get_x()
                y = grids[str(bl)].get_y()
                ax.plot(x,y,color='green',linewidth=0.5)
                
    return grids

# # New Decomposition

# #print('hiii')
# up = 0
# right = 1
# down = 2
# left = 3
# convex = 0
# concave = 1
# graph = 0
# edge = 1


# def get_plot(array,type_,title=None):
#     plt.ion()
#     fig,ax = plt.subplots()
#     if type_ == graph:
#         x = np.empty(0)
#         y = np.empty(0)
#         for i in range(len(array)):
#             x,y = np.append(x,array[i][0]),np.append(y,array[i][1])
#             ax.plot(x[-2:],y[-2:],color='green',linewidth=3)
#             plt.show()
#             plt.pause(0.01)
            
#     elif type_ == edge:
#         count = 0
#         x = np.empty(0)
#         y = np.empty(0)
#         for elem in array:
#             x,y = np.append(x,elem[0]),np.append(y,elem[1])
#             count+=1
#             if count == 2:
#                 ax.plot(x,y)
#                 x = np.empty(0)
#                 y = np.empty(0)
#                 count = 0
#     ax.set_title(title)
#     plt.show()
#     return ax

# def plot_show():
#     plt.figure(figsize = (5,5))
#     plt.show()
    
# def get_coordinates(vertices,boundry_check,action,lim_x,lim_y,edge_length):

#     cur_x = vertices[-1][0]
#     cur_y = vertices[-1][1]

#     boundry = np.array([0,0,0,0])

#     if(cur_x == lim_x):
#         boundry[right] = 1
#         boundry_check[right] = 1
#     if(cur_x == 0):
#         boundry[left] = 1
#         boundry_check[left] = 1
#     if(cur_y == lim_y):
#         boundry[up] = 1
#         boundry_check[up] = 1
#     if(cur_y == 0):
#         boundry[down] = 1
#         boundry_check[down] = 1

#     arr =  (boundry_check == 1)
#     boundry_condition = True

#     for value in arr:
#         if value == False:
#             boundry_condition = False

#     if action == up:
#         next_x = cur_x 
#         next_y = cur_y+1
#         collision = False
#         for i in range(len(vertices)-1):
#             if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
#                 collision=True

#         if(boundry[up] == True):
#             return False,edge_length,vertices

#         elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             edge_length+=1

#             return True,edge_length,vertices

#         elif collision:
#             return False,edge_length,vertices

#         else:
#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             updated_edge_len = edge_length+1

#             trial = 0
#             flag = True
#             actions = [0,1,2,3]
#             while flag:
#                 trial +=1
#                 actions = [0,1,2,3]
#                 choose_action = sample(actions,1)
#                 actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
#                 fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

#                 if fact_vect_temp == True:
#                     indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
#                     unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
#                     if np.array_equal(unique_vertices,temp_vertices[1:]):
#                         return True,edge_len_temp,temp_vertices
#                 if trial == 4:
#                     return False,edge_length,vertices
#                 else:
#                     flag = True

#     elif action == right:
#         next_x = cur_x+1 
#         next_y = cur_y

#         collision = False
#         for i in range(len(vertices)-1):
#             if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
#                 collision=True

#         if(boundry[right] == True):

#             return False,edge_length,vertices

#         elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10  ):

#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             edge_length+=1

#             return True,edge_length,vertices

#         elif(collision):
#             return False,edge_length,vertices

#         else:
#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             updated_edge_len = edge_length+1

#             trial = 0
#             flag = True
#             actions = [0,1,2,3]
#             while flag:
#                 trial +=1
#                 actions = [0,1,2,3]
#                 choose_action = sample(actions,1)
#                 actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
#                 fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

#                 if fact_vect_temp == True:
#                     indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
#                     unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
#                     if np.array_equal(unique_vertices,temp_vertices[1:]):
#                         return True,edge_len_temp,temp_vertices
#                 if trial == 4:
#                     return False,edge_length,vertices
#                 else:
#                     flag = True

#     elif action == down:
#         next_x = cur_x 
#         next_y = cur_y-1
#         collision = False
#         for i in range(len(vertices)-1):
#             if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
#                 collision=True

#         if(boundry[down] == True):
#             return False,edge_length,vertices

#         elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             edge_length+=1

#             return True,edge_length,vertices

#         elif(collision):
#             return False,edge_length,vertices

#         else:
#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             updated_edge_len=edge_length+1                  

#             trial = 0
#             flag = True
#             actions = [0,1,2,3]
#             while flag:
#                 trial +=1
                
#                 choose_action = sample(actions,1)
#                 actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
#                 fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

#                 if fact_vect_temp == True:
#                     indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
#                     unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
#                     if np.array_equal(unique_vertices,temp_vertices[1:]):
#                         return True,edge_len_temp,temp_vertices
#                 if trial == 4:
#                     return False,edge_length,vertices
#                 else:
#                     flag = True

#     elif action == left:
#         next_x = cur_x-1 
#         next_y = cur_y
#         collision = False
#         for i in range(len(vertices)-1):
#             if vertices[i+1][0] == next_x and vertices[i+1][1] == next_y:
#                 collision=True

#         if(boundry[left] == True):
#             return False,edge_length,vertices

#         elif(next_x == vertices[0][0] and next_y == vertices[0][1] and edge_length >= 10 ):

#             vertices = np.append(vertices,[[next_x,next_y]],axis=0)
#             edge_length+=1

#             return True,edge_length,vertices

#         elif(collision):
#             return False,edge_length,vertices

#         else:
#             vertices_ = np.append(vertices,[[next_x,next_y]],axis=0)
#             updated_edge_len=edge_length+1                 

#             trial = 0
#             flag = True
#             actions = [0,1,2,3]
#             while flag:
#                 trial +=1
                
#                 choose_action = sample(actions,1)
#                 actions = list(np.delete(actions, np.where(np.array(actions) == choose_action[0])))
#                 fact_vect_temp,edge_len_temp,temp_vertices = get_coordinates(vertices_,boundry_check, choose_action[0], lim_x,lim_y,updated_edge_len)

#                 if fact_vect_temp == True:
#                     indexes = np.unique(temp_vertices[1:],axis=0, return_index=True)[1]
#                     unique_vertices = [temp_vertices[1:][index] for index in sorted(indexes)]
#                     if np.array_equal(unique_vertices,temp_vertices[1:]):
#                         return True,edge_len_temp,temp_vertices
#                 if trial == 4:
#                     return False,edge_length,vertices
#                 else:
#                     flag = True

# def gen_rectilinear_polygon(lim_x,lim_y):
#     if not type(lim_x) is int:
#         raise TypeError("X should be integer")
#     elif not type(lim_y) is int:
#         raise TypeError("Y should be integer")
#     elif (lim_x) <= 0 :
#         raise Exception("X should be positive integer")
#     elif (lim_y) <= 0 :
#         raise Exception("Y should be positive integer")
#     else:

#         boundry = np.array([0,0,0,0])   #up,right,down,left


#         first_x = randint(0,lim_x-1)
#         first_y = 0

#         vertices = np.array([[first_x,first_y]])

#         junk_1,junk_2,coordinates = get_coordinates(vertices,boundry,right,lim_x,lim_y,0)

#         get_plot(coordinates,graph)

#         return coordinates

# class decomposition:
    
#     def __init__(self):
#         self.main_arr = [[0,0]]
    
#     def get_decomposed_graph(self,array):
        
#         self.main_arr = cp.copy(array)
#         #print(array)
        
#         array = self.get_corner_points(array)
#         #print(array)
        
#         #get_plot(array)

#         H,V = self.get_bipartite(array[::],self.main_arr)

#         #H,V = self.get_absolute_bipartite(H,V)
        
#         remaining_edges = self.get_remaining_edge(array,self.main_arr,H,V)
#         # print(H,V,remaining_edges)
        
#         get_plot(array,graph)
#         get_plot(H,edge)
#         get_plot(V,edge)
#         get_plot(remaining_edges,edge)
#         plot_show()
#     def get_corner_points(self,array):
        
#         x = np.empty(0)
#         y= np.empty(0)
#         for i in range(len(array)):
#             x,y = np.append(x,array[i][0]),np.append(y,array[i][1])

#         deletion_count = 0
#         for i in range(len(x) - 2):
#             j = i - deletion_count

#             pre_x = x[i]
#             recent_x = x[i+1]
#             xi = x[i+2]
#             pre_y = y[i]
#             recent_y = y[i+1]
#             yi = y[i+2]
#             flag = True
#             if pre_x == recent_x == xi:
#                 flag = False
#                 array = np.delete(array,j+1,axis=0)
#                 deletion_count += 1
#             if (pre_y == recent_y == yi) and flag:
#                 flag = True
#                 array = np.delete(array,j+1,axis=0)
#                 deletion_count += 1
                
#         return array
    
#     def get_bipartite(self,coordinates,main_arr):
#         H_of_G = [[0,0]]
#         V_of_G = [[0,0]] 

#         action = [0,0,0,0]

#         angles = self.get_angle(coordinates)
#         coordinates_less_1 = coordinates[:-1]
#         coordinates_less_1_len = len(coordinates_less_1)

#         for i,point in enumerate(coordinates[:-1]):
#             x1 = coordinates[i+1][0]
#             y1 = coordinates[i+1][1]

#             x_dif = x1 - coordinates[i][0]
#             y_dif = y1 - coordinates[i][1]
#             arr = np.append(coordinates_less_1[i:],coordinates_less_1[:-(coordinates_less_1_len-(i))],axis=0)
#             arr1 = arr[1:]
#             arr2 = np.append([arr[-1]],arr[:-1],axis=0)
#             arr2 = arr2[1:]
#             angle_ = np.append(angles[i:],angles[:-(coordinates_less_1_len-(i))],axis=0)
#             angle_1 = angle_[1:]
#             angle_2 = np.append(angle_[-1],angle_[:-1])
#             angle_2 = angle_2[1:]

#             if y_dif >= 1:
#                 action[up] = 1
#                 temp1 = self.get_edge_bipartite(arr1,up,angle_1,main_arr)
#                 temp2 = self.get_edge_bipartite(arr2,down,angle_2,main_arr)
#                 if temp1 != None:
#                     V_of_G = np.append(V_of_G,temp1,axis=0)
#                 elif temp2 != None:
#                     V_of_G = np.append(V_of_G,temp2,axis=0)
#             elif y_dif <= -1:
#                 action[down] = 1
#                 temp1 = self.get_edge_bipartite(arr1,down,angle_1,main_arr)
#                 temp2 = self.get_edge_bipartite(arr2,up,angle_2,main_arr)
#                 if temp1 != None:
#                     V_of_G = np.append(V_of_G,temp1,axis=0)
#                 elif temp2 != None:
#                     V_of_G = np.append(V_of_G,temp2,axis=0)
#             elif x_dif >= 1:
#                 action[right] = 1
#                 temp1 = self.get_edge_bipartite(arr1,right,angle_1,main_arr)
#                 temp2 = self.get_edge_bipartite(arr2,left,angle_2,main_arr)
#                 if temp1 != None:
#                     H_of_G = np.append(H_of_G,temp1,axis=0)
#                 elif temp2 != None:
#                     H_of_G = np.append(H_of_G,temp2,axis=0)
#             elif x_dif <= -1:
#                 action[left] = 1
#                 temp1 = self.get_edge_bipartite(arr1,left,angle_1,main_arr)
#                 temp2 = self.get_edge_bipartite(arr2,right,angle_2,main_arr)
#                 if temp1 != None:
#                     H_of_G = np.append(H_of_G,temp1,axis=0)
#                 elif temp2 != None:
#                     H_of_G = np.append(H_of_G,temp2,axis=0)
#         H,V = self.unique_rows(H_of_G[1:]),self.unique_rows(V_of_G[1:])

#         H,V = self.get_absolute_bipartite(H,V)

#         return H,V
    
        
#     def get_edge_bipartite(self,array,action,angle,main_arr):
#         x,y = array[0][0],array[0][1]
#         if action == up:
#             for i,elem in enumerate(array[1:]):
#                 if elem[0] == x and elem[1] > y:
#                     if angle[i+1] == 1:
#                         if not self.get_obstacle(np.array(array),elem[0],elem[1],up,main_arr):
#                             x1 = x
#                             y1 = y
#                             x2 = elem[0]
#                             y2 = elem[1]
#                             return [[x1,y1],[x2,y2]]
#         elif action == right:
#             for i,elem in enumerate(array[1:]):
#                 if elem[1] == y and elem[0] > x:
#                     if angle[i+1] == 1:
#                         if not self.get_obstacle(np.array(array),elem[0],elem[1],right,main_arr):
#                             x1 = x
#                             y1 = y
#                             x2 = elem[0]
#                             y2 = elem[1]
#                             return [[x1,y1],[x2,y2]]
#         elif action == down:
#             for i,elem in enumerate(array[1:]):
#                 if elem[0] == x and elem[1] < y:
#                     if angle[i+1] == 1:
#                         if not self.get_obstacle(np.array(array),elem[0],elem[1],down,main_arr):
#                             x1 = x
#                             y1 = y
#                             x2 = elem[0]
#                             y2 = elem[1]
#                             return [[x1,y1],[x2,y2]]
#         elif action == left:
#             for i,elem in enumerate(array[1:]):
#                 if elem[1] == y and elem[0] < x:
#                     if angle[i+1] == 1:
#                         if not self.get_obstacle(np.array(array),elem[0],elem[1],left,main_arr):
#                             x1 = x
#                             y1 = y
#                             x2 = elem[0]
#                             y2 = elem[1]
#                             return [[x1,y1],[x2,y2]]


#         return None
    
#     def get_obstacle(self,array,x2,y2,action,main_arr):
#         obstacle = 0
#         x = array[0][0]
#         y = array[0][1]
        
#         if action == up or action == down:
#             #y will be var
#             low = min(y,y2)
#             high = max(y,y2)
#             for j in range(low+1,high):
#                 for i in main_arr[1:]:
#                     if np.array_equal(np.array([x,j]),np.array(i)):
#                         obstacle = 1
                        
#                         break

#         if action == right or action == left:
#             low = min(x,x2)
#             high = max(x,x2)
#             for j in range(low+1,high):
#                 for i in main_arr[1:]:
#                     if np.array_equal(np.array([j,y]),np.array(i)):
#                         obstacle = 1
                        
#                         break
#         return obstacle               
#     #             if i[0] < x and i[1]<=max(y,y2) and i[1]>=min(y,y2):
#     #                 left = 1
#     #             if i[0] > x and i[1]<=max(y,y2) and i[1]>=min(y,y2):
#     #                 right = 1
#     #             if left and right:
#     #                obstacle = 1
#     #                   break

#     #             if i[1] < y and i[0]<=max(x,x2) and i[0]>=min(x,x2):
#     #                 down = 1
#     #             if i[1] > y and i[0]<=max(x,x2) and i[0]>=min(x,x2):
#     #                 up = 1
#     #             if up and down:
#     #                 obstacle = 1
#     #                 break
    
#     def get_absolute_bipartite(self,H,V):
#         V_return = [[0,0]]
#         for i in range(0,len(V),2):
#             high_V = max(V[i][1],V[i+1][1])
#             low_V = min(V[i][1],V[i+1][1])
#             V_x = V[i][0]
#             flag = 0
#             for j in range(0,len(H),2):
#                 high_H = max(H[j][0],H[j+1][0])
#                 low_H = min(H[j][0],H[j+1][0])
#                 H_y = H[j][1]

#                 if H_y < high_V and H_y > low_V:

#                     if V_x < high_H and V_x > low_H:

#                         flag = 1
#             if flag == 0:
#                 V_return = (np.append(V_return,[V[i],V[i+1]],axis=0)).tolist()
#         return H,list(V_return[1:])
    
#     def get_remaining_edge(self,array,main_arr,H,V):
#         main_H = [[0,0]]
#         main_V = [[0,0]]
#         for i in range(0,len(H),2):
#             low = min(H[i][0],H[i+1][0])
#             high = max(H[i][0],H[i+1][0])
            
#             for j in range(low,high+1):
#                 main_H = list(np.append(main_H,[[j,H[i][1]]],axis=0))
#         main_H = main_H[1:]
#         for i in range(0,len(V),2):
#             low = min(V[i][1],V[i+1][1])
#             high = max(V[i][1],V[i+1][1])
            
#             for j in range(low,high+1):
#                 main_V = list(np.append(main_V,[[V[i][0],j]],axis=0))
#         main_V = main_V[1:]
        
#         all_points = main_arr
#         if len(main_H) > 0:
#             all_points = list(np.append(all_points,main_H,axis=0))
#         if len(main_V) > 0:
#             all_points = list(np.append(all_points,main_V,axis=0))
        
#         all_points = np.unique(all_points,axis=0)
        
#         angles = self.get_angle(array)
        
#         remaining_point = [[0,0]]
#         for i,elem in enumerate(array[:-1]):
#             edge_flag = 0
#             if angles[i] == 1:
                
#                 for point in H:
#                     if np.array_equal(point,elem):
#                         edge_flag = 1
#                         break
#                 if edge_flag == 0:
                    
#                     for point in V:
#                         if np.array_equal(point,elem):
#                             edge_flag = 1
#                             break
#                 if edge_flag == 0:
#                     remaining_point = list(np.append(remaining_point,[elem],axis=0))
                    
#         remaining_point = remaining_point[1:]
#         extra_edge = [[0,0]]

#         for point in remaining_point:
# #             print(remaining_point)
# #             print(extra_edge)
#             action = -1
#             for i,elem in enumerate(main_arr):

#                 if np.array_equal(point,elem):
#                     pre = main_arr[i-1]
#                     post = main_arr[i+1]
#                     if pre[1] < elem[1]:
#                         action = up
#                         break
#                     elif pre[1] > elem[1]:
#                         action = down
#                         break
#                     elif pre[1] == elem[1]:
#                         if post[1] > elem[1]:
#                             action = down
#                         else:
#                             action = up

#             if action == up:
#                 extra_edge = np.append(extra_edge,[point],axis=0)
#                 yi = point[1]
#                 xi = point[0]
#                 exit = 0
#                 while(1):
#                     yi+=1
#                     #print(yi)
#                     for j in all_points:
#                         if np.array_equal(j,[xi,yi]):

#                             extra_edge = np.append(extra_edge,[j],axis=0)
#                             exit = 1
#                             break
#                         else:
#                             pass
#                     if exit == 1:
#                         break

#             if action == down:
#                 extra_edge = np.append(extra_edge,[point],axis=0)
#                 yi = point[1]
#                 xi = point[0]
#                 exit = 0
#                 while(1):
#                     yi-=1
#                     #print(xi,yi)
#                     for j in all_points:
#                         if np.array_equal(j,[xi,yi]):

#                             extra_edge = np.append(extra_edge,[j],axis=0)
#                             exit = 1
#                             break
#                         else:
#                             pass
#                     if exit == 1:
#                         break

#         extra_edge = extra_edge[1:]
#         return list(extra_edge)






#     def unique_rows(self,a):
#         #print(a)
#         for i in range(len(a)):
#             for j in range(i):
#                 if np.array_equal(a[j],a[i]):

#                     a[i]=[-1,-1]
#                     #print(a)
#                     break
#         return_a = [elem.tolist() for i,elem in enumerate(a) if not np.array_equal([-1,-1],elem)]
#         #print('final_a',return_a)
#         return return_a

    
#     def get_angle(self,array):

#         past_action = None
#         post_action = None
#         angle = np.zeros(len(array)-1)
#         action = np.array([0,0])

#         for i,point in enumerate(array[:-1]):

#             if i == 0:
#                 past_point = array[-2]
#                 cur_point = array[0]
#                 post_point = array[1]
#             else:
#                 past_point = array[i-1]
#                 cur_point = array[i]
#                 post_point = array[i+1] 

#             x_dif = -(past_point[0] - cur_point[0])
#             y_dif = -(past_point[1] - cur_point[1])
#             if x_dif >= 1:
#                 action[0] = right
#             elif x_dif <= -1:
#                 action[0] = left
#             elif y_dif >= 1:
#                 action[0] = up
#             elif y_dif <= -1:
#                 action[0] = down

#             x_dif = -(cur_point[0] - post_point[0])
#             y_dif = -(cur_point[1] - post_point[1])
#             if x_dif >= 1:
#                 action[1] = right
#             elif x_dif <= -1:
#                 action[1] = left
#             elif y_dif >= 1:
#                 action[1] = up
#             elif y_dif <= -1:
#                 action[1] = down

#             if (action[0]==right and action[1]==down) or (action[0]==left and action[1]==up) or (action[0]==down and action[1]==left) or (action[0]==up and action[1]==right):
#                 angle[i] = concave
#             else:
#                 angle[i] = convex

#         return angle


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
        # plt.show()
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
        if ax!= None:
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
                p[g].border_grid = 1
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
        # plt.show()
        
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
        if ax != None:
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

        plt.show()
    
    # for g in p:
    #     ax.plot(p[g].get_x(),p[g].get_y(),color='orange')
        # if p[g].passage_l:
        #     ax.plot([p[g].get_x()[ind] for ind in range(-2,1)],[p[g].get_y()[ind] for ind in range(-2,1)],color='white')
        # if p[g].passage_r:
        #     ax.plot(p[g].get_x()[1:3],p[g].get_y()[1:3],color='white')
        # if p[g].passage_t:
        #     ax.plot(p[g].get_x()[2:-1],p[g].get_y()[2:-1],color='white')
        # if p[g].passage_b:
        #     ax.plot(p[g].get_x()[:2],p[g].get_y()[:2],color='white')
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
            centroids.append(grid_graph[str([i,j])].centroid)
            if i-inner_grid_width>=x1:
                grid_graph[str([i,j])].l_nei = str([i-inner_grid_width,j])
            if i+inner_grid_width<x2:
                grid_graph[str([i,j])].r_nei = str([i+inner_grid_width,j])
            if j-inner_grid_height>=y1:
                grid_graph[str([i,j])].b_nei = str([i,j-inner_grid_height])
            if j+inner_grid_height<y2:
                grid_graph[str([i,j])].t_nei = str([i,j+inner_grid_height])
            grid_graph[str([i, j])].rectangle_num = rec_num
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


    rectangular_hilbert_curves = get_space_fill_graph(centroids[0][0],centroids[0][1],centroids[-1][0],centroids[-1][1],inner_grid_height,inner_grid_width,ax)

    for i in range(len(rectangular_hilbert_curves[0])):
        if i-1>=0:
            grid_graph[str([int(rectangular_hilbert_curves[0][i] - inner_grid_width/2), int(rectangular_hilbert_curves[1][i] - inner_grid_height/2)])].path_pre = str([int(rectangular_hilbert_curves[0][i-1] - inner_grid_width/2), int(rectangular_hilbert_curves[1][i-1] - inner_grid_height/2)])
        if i+1 < len(rectangular_hilbert_curves[0]):
            grid_graph[str([int(rectangular_hilbert_curves[0][i] - inner_grid_width/2), int(rectangular_hilbert_curves[1][i] - inner_grid_height/2)])].path_post = str([int(rectangular_hilbert_curves[0][i+1] - inner_grid_width/2), int(rectangular_hilbert_curves[1][i+1] - inner_grid_height/2)])

    return grid_graph,centroids,rectangular_hilbert_curves

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

def get_space_fill_graph(x1,y1,x2,y2,inner_grid_height,inner_grid_width,ax):
    height = (abs(y2-y1)//inner_grid_height) + 1
    width = (abs(x2-x1)//inner_grid_width) + 1

    # print("height : ",height)
    # print("weight : ",width)

    x_arr = [x*inner_grid_width for x, y in gilbert2d(width, height)]
    y_arr = [y*inner_grid_height for x, y in gilbert2d(width, height)]

    x,y = x_arr + x1,y_arr + y1 
    if ax!= None:
        ax.plot(x,y,linewidth = 1,alpha=0.3)

    return [np.array(x),np.array(y)]

class Robot:
    def __init__(self,id,r_type):
        self.id = id
        self.type = r_type
        self.present_loc = None
        self.past_loc = None
        self.next_loc = None
        self.body = None
        self.mode = 'b'
        self.trajectory = []
        self.plott = None
    
    def update_trajectory(self,loc):
        self.trajectory.pop(0)
        self.trajectory = self.trajectory + [loc]

def intruder_(grids,ax):
    intruder = Robot(-1,'i')
    intruder.present_loc = list(grids.items())[np.random.randint(len(grids))][0]
    if ax != None:
        intruder.body = ax.scatter([grids[intruder.present_loc].centroid[0]],[grids[intruder.present_loc].centroid[1]],color='red',s=25)
        ax.add_artist(intruder.body)
    return intruder

# Multi searcher
def searchers(num_robots,grids,num_grids_per_rectangle,hsc,grid_width,grid_height,ax):
    robots = [] #   Spawn Robots

    for g in grids:
        if grids[g].passage:
            robots.append(Robot(len(robots)+1,'g'))
            robots[len(robots)-1].present_loc = g
            # grids[robots[len(robots)-1].present_loc].robot_home = True
            if ax != None:
                robots[len(robots)-1].body = ax.scatter([grids[robots[len(robots)-1].present_loc].centroid[0]],[grids[robots[len(robots)-1].present_loc].centroid[1]],color='slateblue',s=10,alpha=0.4)
                ax.add_artist(robots[len(robots)-1].body)
            # plt.show()

    num_robots -= len(robots)
    for r in range(len(num_grids_per_rectangle)):
        start = int(np.sum(num_grids_per_rectangle[:r]))
        num_to_be_added = max(np.round(num_robots*num_grids_per_rectangle[r]/len(grids)),1)
        end = start + num_to_be_added
        counter = 0
        for i in range(int(start),int(end)):
            robots.append(Robot(len(robots)+1,'s'))
            loc = int(counter*num_grids_per_rectangle[r]/(end-start))
            robots[len(robots)-1].present_loc = str([int(hsc[r][0,loc]-grid_width/2),int(hsc[r][1,loc]-grid_height/2)])

            grids[robots[len(robots)-1].present_loc].robot_home = True
            if ax != None:
                robots[len(robots)-1].body = ax.scatter([grids[robots[len(robots)-1].present_loc].centroid[0]],[grids[robots[len(robots)-1].present_loc].centroid[1]],color='green',s=10)
                ax.add_artist(robots[len(robots)-1].body)
            counter += 1
    return robots

# # Single searcher 
# def searchers(num_robots,grids,num_grids_per_rectangle,hsc,grid_width,grid_height,ax):
#     robots = [] #   Spawn Robots

#     start = int(np.sum(num_grids_per_rectangle[:r]))
#     num_to_be_added = max(np.round(num_robots*num_grids_per_rectangle[r]/len(grids)),1)
#     end = start + num_to_be_added
#     counter = 0

#     robots.append(Robot(len(robots)+1,'s'))
#     loc = int(counter*num_grids_per_rectangle[r]/(end-start))
#     robots[len(robots)-1].present_loc = str([int(hsc[r][0,loc]-grid_width/2),int(hsc[r][1,loc]-grid_height/2)])

#     grids[robots[len(robots)-1].present_loc].robot_home = True
#     robots[len(robots)-1].body = ax.scatter([grids[robots[len(robots)-1].present_loc].centroid[0]],[grids[robots[len(robots)-1].present_loc].centroid[1]],color='green',s=2)
#     ax.add_artist(robots[len(robots)-1].body)
#     counter += 1
#     return robots

# k searcher 
# def searchers(num_robots,grids,num_grids_per_rectangle,hsc,grid_width,grid_height,ax):
#     robots = [] #   Spawn Robots

#     for r in range(len(num_grids_per_rectangle)):
#         start = int(np.sum(num_grids_per_rectangle[:r]))
#         num_to_be_added = max(np.round(num_robots*num_grids_per_rectangle[r]/len(grids)),1)
#         end = start + num_to_be_added
#         counter = 0
#         for i in range(int(start),int(end)):
#             robots.append(Robot(len(robots)+1,'s'))
#             loc = int(counter*num_grids_per_rectangle[r]/(end-start))
#             robots[len(robots)-1].present_loc = str([int(hsc[r][0,loc]-grid_width/2),int(hsc[r][1,loc]-grid_height/2)])

#             grids[robots[len(robots)-1].present_loc].robot_home = True
#             if ax is not None:
#                 robots[len(robots)-1].body = ax.scatter([grids[robots[len(robots)-1].present_loc].centroid[0]],[grids[robots[len(robots)-1].present_loc].centroid[1]],color='green',s=10)
#                 ax.add_artist(robots[len(robots)-1].body)
#             counter += 1
#     return robots


# Multi searcher 
static_intruder = 0
if static_intruder:
    path = os.getcwd()
    performance = []
    fig,ax = plt.subplots()
    plt.box(False)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    itera = np.random.randint(50)
    print('itera :',itera)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(20,g_size=10,ax=ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    # plt.pause(1)
    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        # ax = None
        grid_graph,centroids, hsc = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)

        k_max += len(centroids)
        area_reactangles.append(len(centroids))
        HSC = HSC + [np.array(hsc)]
    # plt.ioff()
    # plt.show()
    # plt.pause(1)
    data = []
    
    min_robo_req = 0
    for g in grid_graph:
        if grid_graph[g].passage:
            min_robo_req += 1
    
    k_min = len(area_reactangles)# + min_robo_req
    for k in range(k_min,k_max+1):
        print(k)
        
        
        searchers_time_num = {}#    key = searcher number
        data_run = []
        for runs in range(100):
            # print('num_robos',k)
            for g in grid_graph:
                if grid_graph[g].robot_home:
                    grid_graph[g].robot_home = False
            intruder = intruder_(grid_graph,ax)
            robots = searchers(k,grid_graph,area_reactangles,HSC,grid_width,grid_height,ax)
            intruder_found_state = False
            t = 0
            while not intruder_found_state:
                for rb in robots:
                    if rb.present_loc == intruder.present_loc:
                        intruder_found_state = True
                    if rb.type == 's':
                        rb.past_loc = rb.present_loc
                        if grid_graph[rb.present_loc].robot_home == True:# or grid_graph[rb.present_loc].path_post == None or grid_graph[rb.present_loc].path_pre == None:
                            if rb.mode == 'f':
                                rb.mode = 'b'
                            elif rb.mode == 'b':
                                rb.mode = 'f'
                        elif rb.mode == 'f' and grid_graph[rb.present_loc].path_post == None:
                            rb.mode = 'b'
                        elif rb.mode == 'b' and grid_graph[rb.present_loc].path_pre == None:
                            rb.mode = 'f'
                        if rb.mode == 'f':
                            if grid_graph[rb.present_loc].path_post != None:
                                rb.present_loc = grid_graph[rb.present_loc].path_post
                        else:
                            if grid_graph[rb.present_loc].path_pre!= None:
                                rb.present_loc = grid_graph[rb.present_loc].path_pre
                        rb.body.set_offsets([[grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]]])
                        if len(rb.trajectory)<5:
                            rb.trajectory.append([grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]])
                            if len(rb.trajectory)==4:
                                rb.plott = ax.plot(np.array(rb.trajectory)[:,0],np.array(rb.trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))
                        else:
                            rb.update_trajectory([grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]])
                            rb.plott[0].set_data(np.array(rb.trajectory)[:,0],np.array(rb.trajectory)[:,1])
                ax.set_xlim(40)
                ax.set_ylim(40)
                plt.show()
                plt.pause(0.000001)
                t += 1
                if intruder_found_state:
                    intruder.body.set_visible(False)
                    for rb in robots:
                        rb.body.set_visible(False)
                        if rb.plott != None:
                            rb.plott[0].set_visible(False)
                    searchers_time_num[len(robots)] = t - 1
                    break
        #     data_run.append([k,t])
        # data.append(data_run)
        for rb in robots:
            rb.body.set_visible(False)
            if rb.plott != None:
                rb.plott[0].set_visible(False)
    # fileObject = open(path+'/results/MRIS_s_sfc', 'wb')
    # pkl.dump(data,fileObject)
    # fileObject.close()

dynamic_intruder = 0
if dynamic_intruder:
    path = os.getcwd()
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
    performance = []
    fig,ax = plt.subplots()
    plt.box(False)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    itera = np.random.randint(50)
    print('itera :',itera)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(20,g_size=10,ax=ax)

    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=None)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())

    # plt.pause(1)
    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        # ax = None
        grid_graph,centroids, hsc = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)

        k_max += len(centroids)
        area_reactangles.append(len(centroids))
        HSC = HSC + [np.array(hsc)]
    # plt.ioff()
    # plt.show()
    # plt.pause(1)
    data = []
    min_robo_req = 0
    for g in grid_graph:
        if grid_graph[g].passage:
            min_robo_req += 1
    k_min = len(area_reactangles)# + min_robo_req
    for k in range(k_min,k_max+1):
        print(k)

        searchers_time_num = {}#    key = searcher number
        data_run = []
        for runs in range(100):
            for g in grid_graph:
                if grid_graph[g].robot_home:
                    grid_graph[g].robot_home = False
            intruder = intruder_(grid_graph,ax)
            robots = searchers(k,grid_graph,area_reactangles,HSC,grid_width,grid_height,ax)
            intruder_found_state = False
            t = 0
            while not intruder_found_state:
                intruder_next = np.random.choice(sides)
                intruder.past_loc = intruder.present_loc
                intruder.present_loc = get_neighbour(grid_graph[intruder.present_loc],intruder_next)
                if intruder.present_loc == None:
                    intruder.present_loc = intruder.past_loc 
                intruder.body.set_offsets([[grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]]])
                if len(intruder.trajectory)<10:
                    intruder.trajectory.append([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                    if len(intruder.trajectory)==9:
                        intruder.plott = ax.plot(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1],color='red',linewidth=3,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))
                else:
                    intruder.update_trajectory([grid_graph[intruder.present_loc].centroid[0],grid_graph[intruder.present_loc].centroid[1]])
                    intruder.plott[0].set_data(np.array(intruder.trajectory)[:,0],np.array(intruder.trajectory)[:,1])

                for rb in robots:
                    if rb.present_loc == intruder.present_loc:
                        intruder_found_state = True
                    if rb.type == 's':
                        rb.past_loc = rb.present_loc
                        if grid_graph[rb.present_loc].robot_home == True:# or grid_graph[rb.present_loc].path_post == None or grid_graph[rb.present_loc].path_pre == None:
                            if rb.mode == 'f':
                                rb.mode = 'b'
                            elif rb.mode == 'b':
                                rb.mode = 'f'
                        elif rb.mode == 'f' and grid_graph[rb.present_loc].path_post == None:
                            rb.mode = 'b'
                        elif rb.mode == 'b' and grid_graph[rb.present_loc].path_pre == None:
                            rb.mode = 'f'
                        if rb.mode == 'f':
                            if grid_graph[rb.present_loc].path_post != None:
                                rb.present_loc = grid_graph[rb.present_loc].path_post
                        else:
                            if grid_graph[rb.present_loc].path_pre!= None:
                                rb.present_loc = grid_graph[rb.present_loc].path_pre
                        rb.body.set_offsets([[grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]]])
                        if len(rb.trajectory)<5:
                            rb.trajectory.append([grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]])
                            if len(rb.trajectory)==4:    
                                rb.plott = ax.plot(np.array(rb.trajectory)[:,0],np.array(rb.trajectory)[:,1],color='green',linewidth=2,alpha = 0.3)#,marker='o',markersize=0.1*len(rb.trajectory))

                        else:
                            rb.update_trajectory([grid_graph[rb.present_loc].centroid[0],grid_graph[rb.present_loc].centroid[1]])
                            rb.plott[0].set_data(np.array(rb.trajectory)[:,0],np.array(rb.trajectory)[:,1])
                        
                ax.set_xlim(40)
                ax.set_ylim(40)
                plt.show()
                plt.pause(0.000001)
                t += 1
                if intruder_found_state:
                    intruder.body.set_visible(False)
                    if intruder.plott != None:
                        intruder.plott[0].set_visible(False)
                    for rb in robots:
                        rb.body.set_visible(False)
                        if rb.plott != None:
                            rb.plott[0].set_visible(False)
                    searchers_time_num[len(robots)] = t - 1
                    break
        #     data_run.append([k,t])
        # data.append(data_run)
        intruder.body.set_visible(False)
        if intruder.plott != None:
            intruder.plott[0].set_visible(False)
        for rb in robots:
            rb.body.set_visible(False)
            if rb.plott != None:
                rb.plott[0].set_visible(False)
    # fileObject = open(path+'/results/MRIS_d_sfc', 'wb')
    # pkl.dump(data,fileObject)
    # fileObject.close()

test = 1
if test:
    path = os.path.realpath(os.path.dirname(__file__)) 
    performance = []
    plt.ion()
    fig,ax = plt.subplots()
    plt.box(False)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    itera = np.random.randint(50)
    print('itera :',itera)
    boundary,edges,start_edge,end_edge,grids = Inflate_Cut_algorithm(20,g_size=10,ax=ax)
    
    print(len(grids))
    
    boundary1 = get_area_vertices(np.concatenate((boundary[2:],boundary[:3]),axis=0))

    # D = decomposition()
    # D.get_decomposed_graph(boundary)
    
    # until here boundary contains polygon, boundary1 has polygon vetex set
    # Grid graph formation:
    grids = Grid_graph(boundary1,grid_size=5,ax=ax)
    plt.ioff()
    plt.show()
    rectangles_nodes,grids = get_rectangles(grids,get_area_vertices(boundary),ax=ax)

    rectangles = []
    for r in rectangles_nodes:
        rectangles.append(r.get_corners())
        
    # plt.pause(1)
    
    k_max = 0
    grid_graph = {}
    HSC = []
    grid_height = 5
    grid_width = 5
    area_reactangles = []
    for r in range(len(rectangles)):
        # ax = None
        grid_graph,centroids, hsc = make_grid(min(np.array(rectangles[r])[:,0]),min(np.array(rectangles[r])[:,1]),max(np.array(rectangles[r])[:,0]),max(np.array(rectangles[r])[:,1]),grid_height,grid_width,ax,grid_graph,grids,r)

        k_max += len(centroids)
        area_reactangles.append(len(centroids))
        HSC = HSC + [np.array(hsc)]
    plt.ioff()
    plt.show()
    for i in grid_graph:
        x = grid_graph[i].get_x()
        y = grid_graph[i].get_y()
        ax.plot(x,y,color='orange')
    ax.plot(boundary[:,0],boundary[:,1],color = 'indigo')
    plt.ioff()
    # fig.savefig(path+'/results/5s_grided.pdf',format = "pdf",bbox_inches="tight",pad_inches=0)
    plt.show()

# # Plot the polygons
# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_xlim([-1, 6])
# ax.set_ylim([-1, 6])
# ax.set_xticks(range(-1, 6))
# ax.set_yticks(range(-1, 6))

# # Plot the original polygon
# polygon_patch = plt.Polygon([polygon.vertices[v].coords for v in range(len(polygon.vertices))], fc='none', ec='black')
# ax.add_patch(polygon_patch)

# # Plot the subpolygons
# for subpolygon in subpolygons:
#     subpolygon_patch = plt.Polygon([subpolygon.vertices[v].coords for v in range(len(subpolygon.vertices))], fc='none', ec='red')
#     ax.add_patch(subpolygon_patch)

# # Plot the unused rectangles
# for rectangle in unused_rectangles:
#     rectangle_patch = plt.Rectangle((0, 0), rectangle.width, rectangle.height, fc='none', ec='blue')
#     ax.add_patch(rectangle_patch)

# plt.show()
