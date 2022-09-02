import matplotlib.pyplot as plt
import numpy as np

class Node:
    def __init__(self,id,arena_w,arena_h) -> None:
        self.index = id
        self.node_id = np.array([int(id/arena_w),id%arena_w])
        self.probability = 1/(arena_h*arena_w)
        self.cost = 0
        self.distance = None
        self.parent = None
        self.f = 0
        self.neighbors = []
    
    def add_neighbors(self,grid,column,row):

        neighbor_x = self.node_id[0]
        neighbor_y = self.node_id[1]
    
        if neighbor_x < column - 1:
            self.neighbors.append(grid[neighbor_x+1][neighbor_y])
        if neighbor_x > 0:
            self.neighbors.append(grid[neighbor_x-1][neighbor_y])
        if neighbor_y < row:
            self.neighbors.append(grid[neighbor_x][neighbor_y +1])
        if neighbor_y > 0: 
            self.neighbors.append(grid[neighbor_x][neighbor_y-1])

class Robot:
    def __init__(self,id,arena_w,arena_h,nodes):
        self.id = id
        self.past_loc = None
        self.arena_w = arena_w
        self.arena_h = arena_h
        self.present_loc = self.loc(nodes)
        self.next_loc = None
        self.path_trace = [self.present_loc]
        self.body = None
        self.path_trace_past = None

    def loc(self,nodes):
        probs = probabilities_(nodes,self.arena_w,self.arena_h)
        probs = probs.reshape(self.arena_h*self.arena_w).astype(float)
        random_loc_num = np.random.choice(range(self.arena_h*self.arena_w),p=probs/np.sum(probs))
        return np.array([int(random_loc_num/self.arena_w),random_loc_num%self.arena_w])


def path_finder(origin,goal,grid,fig,ax):
    ax.clear()
    prob_cost_map = np.ones(grid.shape)*0.5
    prob_cost_map[origin[0],origin[1]] = 0.0
    cs = ax.pcolormesh(range(arena_w),range(arena_h),prob_cost_map,shading='auto')
    cs.set_clim(0,1)
    ax.set_aspect('equal', 'box')
    major_ticks = np.arange(0, 101, 1)
    minor_ticks = np.arange(0, 101, 1)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks)
    ax.set_yticks(minor_ticks, minor=True)
    ax.grid(b=True, which='major', color='black', linestyle='-',linewidth = 0.2)
    ax.grid(b=True, which='minor', color='black', linestyle='-',linewidth = 0.2)
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])
    plt.show()
     
    finished = False
    present_loc = origin
    count=0
    probabilities = np.copy(grid)
    path = [present_loc]
    ax.scatter([goal[1]+0.5],[goal[0]+0.5],edgecolors='black', facecolors='none',s=5)
    while not finished:
        next_finder = np.zeros(grid.shape)
        ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],edgecolors='orange', facecolors='none',s=5)
        plt.show()
        left = present_loc - np.array([0,1])
        right = present_loc + np.array([0,1])
        top = present_loc + np.array([1,0])
        bottom = present_loc - np.array([1,0])

        if left[0]>=0 and left[0]<grid.shape[0] and left[1]>=0 and left[1]<grid.shape[1]:
            next_finder[left[0],left[1]] = probabilities[left[0],left[1]].probability
        if right[0]>=0 and right[0]<grid.shape[0] and right[1]>=0 and right[1]<grid.shape[1]:
            next_finder[right[0],right[1]] = probabilities[right[0],right[1]].probability
        if top[0]>=0 and top[0]<grid.shape[0] and top[1]>=0 and top[1]<grid.shape[1]:
            next_finder[top[0],top[1]] = probabilities[top[0],top[1]].probability
        if bottom[0]>=0 and bottom[0]<grid.shape[0] and bottom[1]>=0 and bottom[1]<grid.shape[1]:
            next_finder[bottom[0],bottom[1]] = probabilities[bottom[0],bottom[1]].probability

        loc = np.where(next_finder == np.max(next_finder))
        choosen = np.random.randint(0,len(loc[0]))
        present_loc = np.array([loc[0][choosen],loc[1][choosen]])
        path.append(present_loc)
        
        probabilities[present_loc[0],present_loc[1]].probability = max(0,probabilities[present_loc[0],present_loc[1]].probability - 1/(grid.shape[0]*grid.shape[1]))

        if present_loc[0] == goal[0] and present_loc[1] == goal[1]:
            finished=True
        count += 1

    return path


def to_be_visited(num_grids,visited_and_distance):
    v = -10
    for i in range(num_grids):
        if visited_and_distance[i][0] == 0 and (v < 0 or visited_and_distance[i][1] <= visited_and_distance[v][1]):
            v = i
    return v

def cost_map_(grids):
    cost = np.zeros_like(grids)
    for i in range(grids.shape[0]):
        for j in range(grids.shape[1]):
            cost[i,j] = grids[i,j].cost
    return cost

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
                # ax1.scatter([current.node_id[1]+0.5],[current.node_id[0]+0.5],color='brown',s=50)
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

    print(path)
    # visited_and_distance = [[0,np.Infinity,cost_map[int(i/grids.shape[1]),i%grids.shape[1]]] for i in range(num_grids)]
    # visited_and_distance[source[0]*grids.shape[1]+source[1]] = [0,0,cost_map[int(i/grids.shape[1]),i%grids.shape[1]]]

    # for vertex in range(num_grids):
    #     to_visit = to_be_visited(num_grids,visited_and_distance)
    #     for neighbor in range(num_grids):
    #         if vertices[to_visit][neighbor] == 1 and visited_and_distance[neighbor][0] == 0:
    #             new_cost = visited_and_distance[to_visit][1] + edges[to_visit][neighbor] + visited_and_distance[to_visit][2]
    #             if visited_and_distance[neighbor][1] > new_cost:
    #                 visited_and_distance[neighbor][1] = new_cost
    #         visited_and_distance[to_visit][0] = 1
    # print(visited_and_distance)
    return path

def probabilities_(grids,arena_w,arena_h):
    prob = np.zeros_like(grids)
    for i in range(grids.shape[0]):
        for j in range(grids.shape[1]):
            prob[i,j] = grids[i,j].probability
    return prob

def intruder_(arena_w,arena_h,grids,ax):
    intruder = Robot(-1,arena_w,arena_h,grids)    #   Spawn Intruder
    random_loc_num = np.random.choice(range(arena_h*arena_w))
    intruder.present_loc = np.array([int(random_loc_num/arena_w),random_loc_num%arena_w])
    intru = ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],color='red',s=5)
    ax.add_artist(intru)
    return intruder

def searchers(num_robots,arena_w,arena_h,grids,ax):
    robots = [] #   Spawn Robots
    for i in range(num_robots):
        robots.append(Robot(i,arena_w,arena_h,grids))
        robots[i].body = ax.scatter([robots[i].present_loc[1]+0.5],[robots[i].present_loc[0]+0.5],color='violet',s=5)
        ax.add_artist(robots[i].body)
    return robots

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
                ax.scatter([r.path_trace_past[0][1]+0.5],[r.path_trace_past[0][0]+0.5],edgecolors='green', facecolors='none')
                ax.scatter([p[1]+0.5],[p[0]+0.5],edgecolors='yellow', facecolors='none',s=5)
                ax.scatter([r.path_trace_past[-1][1]+0.5],[r.path_trace_past[-1][0]+0.5],edgecolors='black', facecolors='none')

def update_prob_costs(grids,present_locs,arena_h,arena_w):
    for row in range(grids.shape[0]): #   Update the Costs and Probabilities of the grids
        for col in range(grids.shape[1]):
            loc = np.array([row,col])
            for l in present_locs:
                if np.linalg.norm(loc-np.array(l))==0:
                    grids[loc[0],loc[1]].probability -= 1/(arena_h*arena_w)
                    if grids[loc[0],loc[1]].probability<0:
                        grids[loc[0],loc[1]].probability=0
                    grids[loc[0],loc[1]].cost += 1
                else:
                    grids[loc[0],loc[1]].probability += 1*len(present_locs)/(arena_h*arena_w*(arena_h*arena_w))

def plot_mesh(fig,ax,arena_w,arena_h,map,without_bar=0):
    cs = ax.pcolormesh(range(arena_w),range(arena_h),map,shading='auto')
    # neg = ax.imshow(map, cmap='Reds_r', interpolation='none')
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
    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

arena_w = 10
arena_h = 10

grids = np.array([Node(i,arena_h=arena_h,arena_w=arena_w) for i in range(arena_w*arena_h)]).reshape((arena_h,arena_w))
probabilities = probabilities_(grids,arena_w,arena_h)

plt.ion()
fig,ax = plt.subplots()

plot_mesh(fig,ax,arena_w,arena_h,probabilities)

time = []
agents = []

for num_robots in range(5,25,5):    
    t = 0   #   Initialize time
    search_state =  0 #   0 = not found, 1 = found
    
    robots = searchers(num_robots,arena_w,arena_h,grids,ax)
    intruder = intruder_(arena_w,arena_h,grids,ax)
    ax.clear()
    plot_mesh(fig,ax,arena_w,arena_h,probabilities,without_bar = 1)
    while not search_state: #   While Intruder is not found
        # ax.clear()
        present_locs = visited_loc(robots)

        update_prob_costs(grids,present_locs,arena_h,arena_w)

        probabilities = probabilities_(grids,arena_w,arena_h)

        intru = ax.scatter([intruder.present_loc[1]+0.5],[intruder.present_loc[0]+0.5],color='red',s=5)
        ax.add_artist(intru)
        cs = ax.pcolormesh(range(arena_h),range(arena_w),np.array(probabilities),shading='auto')
        ax.set_aspect('equal', 'box')
        major_ticks = np.arange(0, 101, 1)
        minor_ticks = np.arange(0, 101, 1)
        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)
        ax.grid(b=True, which='major', color='black', linestyle='-',linewidth = 0.2)
        ax.grid(b=True, which='minor', color='black', linestyle='-',linewidth = 0.2)
        ax.axes.xaxis.set_ticks([])
        ax.axes.yaxis.set_ticks([])
        for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
            robots[i].body = ax.scatter([robots[i].present_loc[1]+0.5],[robots[i].present_loc[0]+0.5],color='violet',s=5)
            ax.add_artist(robots[i].body)
        plot_past_path(ax,robots)
        plt.show()
        plt.pause(0.001)

        paths = []
        for i in range(num_robots): #   Traverse each robot to new location from its past location and check presence of Intruder along the path
            robots[i].next_loc = robots[i].loc(grids)

            for g in grids:
                for g_ij in g:
                    g_ij.parent = None
            path = Astar(grids,robots[i].present_loc,robots[i].next_loc)#path_finder(robots[i].present_loc,robots[i].next_loc,grids,fig,ax)

            for p in range(len(path)-1):
                ax.scatter([path[p][1]+0.5],[path[p][0]+0.5],color='yellow',s=5)
                ax.plot([path[p][1]+0.5,path[p+1][1]+0.5],[path[p][0]+0.5,path[p+1][0]+0.5],color= 'black',alpha=0.2)
            ax.scatter([path[0][1]+0.5],[path[0][0]+0.5],edgecolors='white', facecolors='none')
            ax.scatter([path[-1][1]+0.5],[path[-1][0]+0.5],marker='x',color = 'black')
            plt.show()
            plt.pause(0.001)
            
            robots[i].path_trace_past = path#robots[i].path_trace
            robots[i].path_trace = path
            paths+= path

            robots[i].past_loc = robots[i].present_loc
            robots[i].present_loc = robots[i].next_loc

        t += 1
        for p in paths:
            if np.linalg.norm(intruder.present_loc-np.array(p))==0:
                search_state = 1
                break       
        
        

    time.append(t)  #   Store time taken to find the intruder
    agents.append(num_robots)   #   Store number of robots utilized for search of static intruder

plt.ioff()
fig,ax = plt.subplots()
ax.plot(time,agents)
ax.set_xlabel("Time taken")
ax.set_ylabel("Number of search agents")
plt.show()