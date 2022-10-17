from audioop import avg
import matplotlib.pyplot as plt
import numpy as np
import faulthandler
import copy as cp

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
        if not isinstance(r.path_trace,type(None)):
            present_locs += r.path_trace
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
                    
                    robots[i].path_progressor += 1
                    robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment

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
                    
                    robots[i].path_progressor += 1
                    robots[i].present_loc = robots[i].path_trace[robots[i].path_progressor] # add index increment

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