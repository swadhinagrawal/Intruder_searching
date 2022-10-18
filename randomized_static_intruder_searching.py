from audioop import avg
import matplotlib.pyplot as plt
import numpy as np
import faulthandler
import copy as cp

# Random simple rectilinear polygon generator

def Inflate_Cut_algorithm(num_vertices,ax=None):
    '''
    Input: Num_vertices (Must be even and >=4)
    Returns: Array of vertices forming the random simple rectilinear region
    '''
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

class Node:
    def __init__(self,id,arena_w,arena_h):
        self.node_id = np.array([int(id/arena_w),id%arena_w])
        self.probability = 1.0/(arena_h*arena_w)
        self.cost = 0.0
        self.distance = None
        self.parent = None
        self.f = 0.0

class Robot:
    def __init__(self,id,arena_w,arena_h,nodes,r_type):
        self.id = id
        self.type = r_type
        self.past_loc = None
        self.arena_w = arena_w
        self.arena_h = arena_h
        self.start = 1
        self.present_loc = self.loc(nodes)
        self.next_loc = None
        self.path_trace = [self.present_loc]
        self.body = None
        self.path_trace_past = None
        self.arena = cp.copy(nodes)
        self.path_progressor = 0
        

    def loc(self,nodes):
        if self.type == 's':
            probs = probabilities_(nodes)
            probs = probs.reshape(self.arena_h*self.arena_w).astype(float)
            random_loc_num = np.random.choice(range(self.arena_h*self.arena_w),p=probs/np.sum(probs))
            return np.array([int(random_loc_num/self.arena_w),random_loc_num%self.arena_w])
        else:
            if self.start:
                random_loc_num = np.random.choice(range(self.arena_h*self.arena_w))
                self.start = 0
                return np.array([int(random_loc_num/self.arena_w),random_loc_num%self.arena_w])
            else:
                multiplier = np.array([[0,1],[1,0]])
                choose_multiplier = np.random.randint(0,2)
                rand_loc = np.random.randint(-1,2,size=2)* multiplier[choose_multiplier]
                new_loc = self.present_loc - rand_loc
                # if new_loc[0]>=self.arena_h or new_loc[1]>=self.arena_w:
                #     new_loc = self.present_loc
                while new_loc[0]>=self.arena_h or new_loc[1]>=self.arena_w or new_loc[0]<=0 or new_loc[1]<=0:
                    choose_multiplier = np.random.randint(0,2)
                    rand_loc = np.random.randint(-1,2,size=2)* multiplier[choose_multiplier]
                    new_loc = self.present_loc - rand_loc
                return new_loc

def cost_map_(grids):
    cost = np.zeros_like(grids)
    for i in range(grids.shape[0]):
        for j in range(grids.shape[1]):
            cost[i,j] = grids[i,j].cost
    return cost

def probabilities_(grids):
    prob = np.zeros_like(grids)
    for i in range(grids.shape[0]):
        for j in range(grids.shape[1]):
            prob[i,j] = grids[i,j].probability
    return prob

def Astar(grids,source,destination):
    # plt.ion()
    # fig1,ax1 = plt.subplots()
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
    edges = np.copy(vertices)
    temp_grids = np.copy(grids)
    cost_map = cost_map_(grids)
    # plot_mesh(fig1,ax1,grids.shape[1],grids.shape[0],cost_map)
    num_grids = len(vertices[0])
    temp_grids[source[0],source[1]].f = 0
    temp_grids[source[0],source[1]].distance = 0
    temp_grids[destination[0],destination[1]].f = 0
    open_nodes = [temp_grids[source[0],source[1]]]
    closed_nodes = []

    while len(open_nodes) > 0:
        current_node = open_nodes[0]
        
        # ax1.scatter([current_node.node_id[1]+0.5],[current_node.node_id[0]+0.5],color = 'green')
        current_index = 0
        for i, item in enumerate(open_nodes):
            if item.f < current_node.f:
                current_node = item
                current_index = i
        # ax1.scatter([current_node.node_id[1]+0.5],[current_node.node_id[0]+0.5],color = 'orange')
        
        # Pop current off open list, add to closed list
        open_nodes.pop(current_index)
        closed_nodes.append(current_node)

        # Found the goal
        if current_node == temp_grids[destination[0],destination[1]]:
            # ax1.scatter([destination[1]+0.5],[destination[0]+0.5],color = 'violet')
            # ax1.scatter([source[1]+0.5],[source[0]+0.5],color = 'black')

            path = []
            current = current_node
            while current is not None:
                path.append(current.node_id)
                # ax1.scatter([current.node_id[1]+0.5],[current.node_id[0]+0.5],color='brown',s=250)
                current = current.parent
                
                
            return path[::-1] # Return reversed path

        # Generate children
        children = []
        for new_position in [(0, -1), (0, 1), (-1, 0), (1, 0)]: # Adjacent squares

            # Get node position
            node_position = (current_node.node_id[0] + new_position[0], current_node.node_id[1] + new_position[1])
            
        
            # Make sure within range
            if (node_position[0] < 0 or node_position[0] > grids.shape[0]-1) or (node_position[1] > grids.shape[1]-1 or node_position[1] < 0):
                continue

            # Create new node
            new_node = temp_grids[node_position[0],node_position[1]]

            # Append
            if np.linalg.norm(np.array(node_position)-source)!=0:
                # ax1.scatter([node_position[1]+0.5],[node_position[0]+0.5],color = 'brown')
                children.append(new_node)

        # Loop through children
        for child in children:
            # ax1.scatter([child.node_id[1]+0.5],[child.node_id[0]+0.5],edgecolors = 'white',facecolor='none')
            # Child is on the closed list
            for closed_child in closed_nodes:
                if child == closed_child:
                    # ax1.scatter([child.node_id[1]+0.5],[child.node_id[0]+0.5],color = 'orange')
                    continue
            
            if child not in closed_nodes:
                # Create the f, g, and h values
                if child==temp_grids[source[0],source[1]]:
                    print('Doomed')
                child.distance = current_node.distance + 1
                child.f = child.distance + child.cost
                child.parent = current_node
                # ax1.plot([child.node_id[1]+0.5,child.parent.node_id[1]+0.5],[child.node_id[0]+0.5,child.parent.node_id[0]+0.5],color= 'black')
                # ax1.scatter([child.parent.node_id[1]+0.5],[child.parent.node_id[0]+0.5],edgecolors= 'red',facecolor='none')

            # Child is already in the open list
            for open_node in open_nodes:
                if child == open_node and child.distance > open_node.distance and child.cost > open_node.cost:
                    # ax1.scatter([child.node_id[1]+0.5],[child.node_id[0]+0.5],color = 'cyan')
                    continue

            # Add the child to the open list
            if child not in open_nodes+closed_nodes:
                open_nodes.append(child)
            # ax1.scatter([destination[1]+0.5],[destination[0]+0.5],color = 'violet')
            # ax1.scatter([source[1]+0.5],[source[0]+0.5],color = 'black')
    return path

def intruder_(arena_w,arena_h,grids,ax):
    intruder = Robot(-1,arena_w,arena_h,grids,'i')    #   Spawn Intruder
    random_loc_num = np.random.choice(range(arena_h*arena_w))
    intruder.present_loc = np.array([int(random_loc_num/arena_w),random_loc_num%arena_w])
    intru = ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],color='red',s=25)
    ax.add_artist(intru)
    return intruder

def intruder_plot(intruder,ax):
    intru = ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],color='red',s=25)
    ax.add_artist(intru)

def searchers(num_robots,arena_w,arena_h,grids,ax):
    robots = [] #   Spawn Robots
    for i in range(num_robots):
        robots.append(Robot(i,arena_w,arena_h,grids,'s'))
        robots[i].body = ax.scatter([robots[i].present_loc[1]+0.5],[robots[i].present_loc[0]+0.5],color='pink',s=25)
        ax.add_artist(robots[i].body)
    return robots

def searchers_plot(robots,ax):
    for r in robots:
        r.body = ax.scatter([r.present_loc[1]+0.5],[r.present_loc[0]+0.5],color='pink',s=25)
        ax.add_artist(r.body)

def visited_loc(robots):
    present_locs = []
    for r in robots:
        present_locs += [list(r.present_loc)]
    return present_locs

def plot_past_path(ax,robots):
    for r in robots:
        if not isinstance(r.path_trace_past,type(None)):
            for p in r.path_trace_past:
                ax.scatter([p[1]+0.5],[p[0]+0.5],edgecolors='maroon', facecolors='none',s=60)
            ax.scatter([r.path_trace_past[0][1]+0.5],[r.path_trace_past[0][0]+0.5],edgecolors='black',s=25, facecolors='none')
            ax.scatter([r.path_trace_past[-1][1]+0.5],[r.path_trace_past[-1][0]+0.5],edgecolors='white',s=25, facecolors='none')

def update_prob_costs(grids,present_locs,arena_h,arena_w):
    for row in range(grids.shape[0]): #   Update the Costs and Probabilities of the grids
        for col in range(grids.shape[1]):
            loc = np.array([row,col])
            for l in present_locs:
                if np.linalg.norm(loc-np.array(l))==0:
                    grids[loc[0],loc[1]].probability -= 1/(arena_h*arena_w)
                    if grids[loc[0],loc[1]].probability<0:
                        grids[loc[0],loc[1]].probability=0
                    grids[loc[0],loc[1]].cost += 1/(arena_h*arena_w)
                else:
                    grids[loc[0],loc[1]].probability += 1/(arena_h*arena_w*len(present_locs))

def plot_mesh(fig,ax,arena_w,arena_h,grid,title,without_bar=0):
    grid = np.round(np.array(grid).astype(float),decimals=4)
    cs = ax.pcolormesh(grid,shading='auto',cmap='plasma')
    if not without_bar:
        cbar = fig.colorbar(cs,orientation='vertical')
        cbar.set_ticks(np.arange(0,1,0.1))
    cs.set_clim(0,1)
    ax.set_aspect('equal', 'box')
    major_ticks = np.arange(0, arena_w+1, 1)
    minor_ticks = np.arange(0, arena_h+1, 1)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    ax.grid(b=True, which='major', color='black', linestyle='-',linewidth = 0.2)
    ax.grid(b=True, which='minor', color='black', linestyle='-',linewidth = 0.2)
    ax.set_title(title)
    # ax.axes.xaxis.set_ticks([])
    # ax.axes.yaxis.set_ticks([])


# Generalize these for generic rectilinear simple polygon
static_intruder = 1
if static_intruder:
    arena_w = 10
    arena_h = 10

    time = []
    agents = []
    plt.ion()
    fig,ax = plt.subplots()
    fig1,ax1 = plt.subplots()
    faulthandler.enable()
    for num_robots in range(1,25,1):
        avg_time = []
        for runs in range(10):
            t = 0   #   Initialize time
            search_state =  0 #   0 = not found, 1 = found
            ax.clear()
            ax1.clear()
            grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))
            probabilities = probabilities_(grids)
            costs = cost_map_(grids)
            if num_robots ==1 and runs ==0:
                without_bar = 0
            else:
                without_bar = 1
            
            plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar)
            plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar)    

            robots = searchers(num_robots,arena_w,arena_h,cp.copy(grids),ax)
            intruder = intruder_(arena_w,arena_h,cp.copy(grids),ax)
            intruder_plot(intruder,ax1)
            searchers_plot(robots,ax1)
            
            while not search_state: #   While Intruder is not found
                ax.clear()
                ax1.clear()
                present_locs = visited_loc(robots)
                update_prob_costs(grids,present_locs,arena_h,arena_w)
                probabilities = probabilities_(grids)
                costs = cost_map_(grids)
                plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar=1)
                plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar=1)

                intruder_plot(intruder,ax1)
                searchers_plot(robots,ax1)

                plot_past_path(ax,robots)
                plt.show()
                plt.pause(0.00001)

                # Intruder update
                intruder.next_loc = intruder.present_loc

                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc

                ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],edgecolors='white',s=25, facecolors='none')
                ax.scatter([intruder.next_loc[1]+0.5],[intruder.next_loc[0]+0.5],marker='x',color = 'green')
                plt.show()
                plt.pause(0.00001)
                

                #   Robots update

                # paths = []
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if np.linalg.norm(robots[i].present_loc-robots[i].path_trace[-1])==0:
                        robots[i].next_loc = robots[i].loc(grids)

                        for g in grids:
                            for g_ij in g:
                                g_ij.parent = None
                        path = Astar(grids,robots[i].present_loc,robots[i].next_loc)
                        for p in range(len(path)-1):
                            ax.scatter([path[p][1]+0.5],[path[p][0]+0.5],edgecolors='cyan',s=25, facecolors='none')

                            # ax.plot([path[p][1]+0.5,path[p+1][1]+0.5],[path[p][0]+0.5,path[p+1][0]+0.5],color= 'black',alpha=0.2)
                        ax.scatter([path[0][1]+0.5],[path[0][0]+0.5],edgecolors='white',s=25, facecolors='none')
                        ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],marker='x',color = 'green')
                        plt.show()
                        plt.pause(0.00001)
                        
                        robots[i].path_trace_past = path#robots[i].path_trace
                        robots[i].path_trace = path
                        robots[i].path_progressor = 0

                    # paths+= list(robots[i].present_loc)

                    robots[i].past_loc = robots[i].present_loc
                    
                    robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment
                    robots[i].path_progressor += 1

                t += 1
                for r in robots:
                    if np.linalg.norm(intruder.present_loc-r.present_loc)==0:
                        search_state = 1
                        break  
                plt.show()
                plt.pause(0.000001)     
            avg_time.append(t)
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder


    plt.ioff()
    fig,ax = plt.subplots()
    ax.plot(time,agents)
    ax.set_xlabel("Time taken")
    ax.set_ylabel("Number of search agents")
    plt.savefig('static_intruder.eps')
    plt.show()

dynamic_intruder = 1
if dynamic_intruder:
    arena_w = 10
    arena_h = 10

    time = []
    agents = []
    plt.ion()
    fig,ax = plt.subplots()
    fig1,ax1 = plt.subplots()
    faulthandler.enable()
    for num_robots in range(1,25,1):
        avg_time = []
        for runs in range(10):
            t = 0   #   Initialize time
            search_state =  0 #   0 = not found, 1 = found
            ax.clear()
            ax1.clear()
            grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))
            probabilities = probabilities_(grids)
            costs = cost_map_(grids)
            if num_robots ==1 and runs ==0:
                without_bar = 0
            else:
                without_bar = 1
            
            plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar)
            plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar)    

            robots = searchers(num_robots,arena_w,arena_h,cp.copy(grids),ax)
            intruder = intruder_(arena_w,arena_h,cp.copy(grids),ax)
            intruder_plot(intruder,ax1)
            searchers_plot(robots,ax1)
            
            while not search_state: #   While Intruder is not found
                ax.clear()
                ax1.clear()
                present_locs = visited_loc(robots)
                update_prob_costs(grids,present_locs,arena_h,arena_w)
                probabilities = probabilities_(grids)
                costs = cost_map_(grids)
                plot_mesh(fig,ax,arena_w,arena_h,probabilities,'Probability Map',without_bar=1)
                plot_mesh(fig1,ax1,arena_w,arena_h,costs,'Cost Map',without_bar=1)

                intruder_plot(intruder,ax1)
                searchers_plot(robots,ax1)

                plot_past_path(ax,robots)
                plt.show()
                plt.pause(0.00001)

                # Intruder update
                intruder.next_loc = intruder.loc(grids)

                intruder.past_loc = intruder.present_loc
                intruder.present_loc = intruder.next_loc

                ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],edgecolors='white',s=25, facecolors='none')
                ax.scatter([intruder.next_loc[1]+0.5],[intruder.next_loc[0]+0.5],marker='x',color = 'green')
                plt.show()
                plt.pause(0.00001)
                

                #   Robots update

                paths = []
                for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
                    if np.linalg.norm(robots[i].present_loc-robots[i].path_trace[-1])==0:
                        robots[i].next_loc = robots[i].loc(grids)

                        for g in grids:
                            for g_ij in g:
                                g_ij.parent = None
                        path = Astar(grids,robots[i].present_loc,robots[i].next_loc)
                        for p in range(len(path)-1):
                            ax.scatter([path[p][1]+0.5],[path[p][0]+0.5],edgecolors='cyan',s=25, facecolors='none')

                            # ax.plot([path[p][1]+0.5,path[p+1][1]+0.5],[path[p][0]+0.5,path[p+1][0]+0.5],color= 'black',alpha=0.2)
                        ax.scatter([path[0][1]+0.5],[path[0][0]+0.5],edgecolors='white',s=25, facecolors='none')
                        ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],marker='x',color = 'green')
                        plt.show()
                        plt.pause(0.00001)
                        
                        robots[i].path_trace_past = path#robots[i].path_trace
                        robots[i].path_trace = path
                        robots[i].path_progressor = 0

                    paths+= list(robots[i].present_loc)

                    robots[i].past_loc = robots[i].present_loc

                    robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment
                    robots[i].path_progressor += 1

                t += 1
                for r in robots:
                    if np.linalg.norm(intruder.present_loc-r.present_loc)==0:
                        search_state = 1
                        break  
                plt.show()
                plt.pause(0.000001)     
            avg_time.append(t)
        time.append(np.sum(avg_time)/len(avg_time))  #   Store time taken to find the intruder
        agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

    plt.ioff()
    fig,ax = plt.subplots()
    ax.plot(time,agents)
    ax.set_xlabel("Time taken")
    ax.set_ylabel("Number of search agents")
    plt.savefig('dynamic_intruder.eps')
    plt.show()


