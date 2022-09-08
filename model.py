import numpy as np
import random as r
import math
from numba import jit, prange
from simulation import Simulation, record_time, template_params
import backend
# from pythonabm import Simulation, record_time, template_params
# from pythonabm.backend import record_time, check_direct, template_params, check_existing, get_end_step, Graph, \
#     progress_bar, starting_params, check_output_dir, assign_bins_jit, get_neighbors_cpu, get_neighbors_gpu


@jit(nopython=True, parallel=True)
def get_neighbor_forces(number_edges, edges, edge_forces, locations, center, types, radius, alpha=10, r_e=1.01,
                        u_bb=5, u_rb=1, u_yb=1, u_rr=20, u_ry=12, u_yy=30, u_repulsion=10000):
    for index in range(number_edges):
        # get indices of cells in edge
        cell_1 = edges[index][0]
        cell_2 = edges[index][1]
        adhesion_values = np.reshape(np.array([u_bb, u_rb, u_yb, u_yb, u_rr, u_ry, u_rb, u_ry, u_yy]), (3, 3))
        # get cell positions
        cell_1_loc = locations[cell_1] - center
        cell_2_loc = locations[cell_2] - center

        # get new location position
        vec = cell_2_loc - cell_1_loc
        dist2 = vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2

        # based on the distance apply force differently
        if dist2 == 0:
            edge_forces[index][0] = 0
            edge_forces[index][1] = 0
        else:
            dist = dist2 ** (1/2)
            if 0 < dist2 < (2 * radius) ** 2:
                edge_forces[index][0] = -1 * u_repulsion * (vec / dist)
                edge_forces[index][1] = 1 * u_repulsion * (vec / dist)
            else:
                # get the cell type
                cell_1_type = types[cell_1]
                cell_2_type = types[cell_2]
                u = adhesion_values[cell_1_type, cell_2_type]
                # get value prior to applying type specific adhesion const
                value = (dist - r_e) * (vec / dist)
                edge_forces[index][0] = u * value
                edge_forces[index][1] = -1 * u * value
    return edge_forces


@jit(nopython=True, parallel=True)
def get_gravity_forces(number_cells, locations, center, well_rad, net_forces, grav=1):
    for index in range(number_cells):
        new_loc = locations[index] - center
        # net_forces[index] = -grav * (new_loc / well_rad) * \
        #                     np.sqrt(1 - (np.linalg.norm(new_loc) / well_rad) ** 2)
        new_loc_sum = new_loc[0] ** 2 + new_loc[1] ** 2 + new_loc[2] ** 2
        net_forces[index] = -grav * (new_loc / well_rad) * new_loc_sum ** (1/2)
    return net_forces


@jit(nopython=True)
def convert_edge_forces(number_edges, edges, edge_forces, neighbor_forces):
    for index in range(number_edges):
        # get indices of cells in edge
        cell_1 = edges[index][0]
        cell_2 = edges[index][1]

        neighbor_forces[cell_1] += edge_forces[index][0]
        neighbor_forces[cell_2] += edge_forces[index][1]

    return neighbor_forces

def seed_cells(num_agents, center, radius):
    theta = 2 * np.pi * np.random.rand(num_agents).reshape(num_agents, 1)
    rad = radius * np.sqrt(np.random.rand(num_agents)).reshape(num_agents, 1)
    x = rad * np.cos(theta) + center[0]
    y = rad * np.sin(theta) + center[1]
    z = np.zeros((num_agents, 1)) + center[2]
    locations = np.hstack((x, y, z))
    return locations

def calculate_rate(combined_percent):
    return 1 - np.power(1-combined_percent, 1/720)

def calculate_transition_rate(combined_percent, t):
    transition_rate = 1 - (1-combined_percent) ** t
    return transition_rate

class TestSimulation(Simulation):
    """ This class inherits the Simulation class allowing it to run a
        simulation with the proper functionality.
    """
    def __init__(self, yaml_file='general.yaml'):
        # initialize the Simulation object
        Simulation.__init__(self)

        # read parameters from YAML file and add them to instance variables
        self.yaml_parameters(yaml_file)
        self.yaml_name = yaml_file
        # aba/dox/cho ratio
        self.cho_ratio = 1 - (self.aba_ratio + self.dox_ratio)
        self.aba_color = np.array([255, 255, 0], dtype=int) #yellow
        self.dox_color = np.array([255, 50, 50], dtype=int) #red
        self.cho_color = np.array([50, 50, 255], dtype=int)


        self.initial_seed_rad = self.well_rad * self.initial_seed_ratio
        self.dim = np.asarray(self.size)
        self.size = self.dim * self.well_rad
        self.center = np.array([self.size[0] / 2, self.size[1] / 2, 0])

    def setup(self):
        """ Overrides the setup() method from the Simulation class.
        """
        # determine the number of agents for each cell type, pre-seeded
        # num_aba = int(self.num_to_start * self.aba_ratio)
        # num_dox = int(self.num_to_start * self.dox_ratio)
        # num_cho = int(self.num_to_start * self.cho_ratio)

        # Seeding cells with stochastic transition rates
        num_aba = 0
        num_dox = 0
        num_cho = int(self.num_to_start)

        # add agents to the simulation
        self.add_agents(num_aba, agent_type="ABA")
        self.add_agents(num_dox, agent_type="DOX")
        self.add_agents(num_cho, agent_type="CHO")

        # indicate agent arrays and create the arrays with initial conditions
        self.indicate_arrays("locations", "radii", "colors", "cell_type", "division_set", "div_thresh")

        # generate random locations for cells
        self.locations = seed_cells(self.number_agents, self.center, self.initial_seed_rad)
        self.radii = self.agent_array(initial=lambda: self.cell_rad)

        # Define cell types, 2 is ABA, 1 is DOX, 0 is non-cadherin expressing cho cells
        self.cell_type = self.agent_array(dtype=int, initial={"ABA": lambda: 2, "DOX": lambda: 1, "CHO": lambda: 0})
        self.colors = self.agent_array(dtype=int, vector=3, initial={"ABA": lambda: self.aba_color, "DOX": lambda: self.dox_color, "CHO": lambda: self.cho_color})

        # setting division times (in seconds):
        # Not used in model
        self.div_thresh = self.agent_array(initial={"ABA": lambda: 1, "DOX": lambda: 1, "CHO": lambda: 1})
        self.division_set = self.agent_array(initial={"ABA": lambda: 0, "DOX": lambda: 0, "CHO": lambda: 0})

        # save parameters to text file
        self.save_params(self.yaml_name)

        #indicate and create graphs for identifying neighbors
        self.indicate_graphs("neighbor_graph")
        self.neighbor_graph = self.agent_graph()

        #stochastic gene circuit
        self.transition_rate = calculate_rate(self.dox_ratio + self.aba_ratio)
        self.transition = np.random.rand(self.number_agents)
        self.cell_fate()

        # record initial values
        self.step_values()
        self.step_image()

    def step(self):
        """ Overrides the step() method from the Simulation class.
        """
        # preform subset force calculations
        for i in range(self.sub_ts):
            # get all neighbors within threshold (1.6 * diameter)
            self.get_neighbors(self.neighbor_graph, self.cell_interaction_rad * self.cell_rad)
            # move the cells and track total repulsion vs adhesion forces
            self.move_parallel()
        # get the following data. We can generate images at each time step, but right now that is not needed.

        # add/remove agents from the simulation
        # self.update_populations()
        print(f'Num_ABA: {len(np.argwhere(self.cell_type == 2))}, Num_dox: {len(np.argwhere(self.cell_type == 1))}')

        # Update cell fate
        self.cell_fate()

        self.step_values()
        self.step_image()
        self.temp()
        self.data()

    def end(self):
        """ Overrides the end() method from the Simulation class.
        """
        self.step_values()
        self.step_image()
        self.create_video()

    @record_time
    def update_populations(self):
        """ Adds/removes agents to/from the simulation by adding/removing
            indices from the cell arrays and any graphs.
        """
        # get indices of hatching/dying agents with Boolean mask
        add_indices = np.arange(self.number_agents)[self.hatching]
        remove_indices = np.arange(self.number_agents)[self.removing]

        # count how many added/removed agents
        num_added = len(add_indices)
        num_removed = len(remove_indices)
        # go through each agent array name
        for name in self.array_names:
            # copy the indices of the agent array data for the hatching agents
            copies = self.__dict__[name][add_indices]

            # add indices to the arrays
            self.__dict__[name] = np.concatenate((self.__dict__[name], copies), axis=0)

            # if locations array
            if name == "locations":
                # go through the number of cells added
                for i in range(num_added):
                    # get mother and daughter indices
                    mother = add_indices[i]
                    daughter = self.number_agents + i

                    # move distance of radius in random direction
                    vec = self.radii[i] * self.random_vector()
                    self.__dict__[name][mother] += vec
                    self.__dict__[name][daughter] -= vec

            # reset division time
            if name == "division_set":
                # go through the number of cells added
                for i in range(num_added):
                    # get mother and daughter indices
                    mother = add_indices[i]
                    daughter = self.number_agents + i

                    # set division counter to zero
                    self.__dict__[name][mother] = 0
                    self.__dict__[name][daughter] = 0

            # set new division threshold
            if name == "division_threshold":
                # go through the number of cells added
                for i in range(num_added):
                    # get daughter index
                    daughter = self.number_agents + i

                    # set division threshold based on cell type
                    self.__dict__[name][daughter] = 1

            # remove indices from the arrays
            self.__dict__[name] = np.delete(self.__dict__[name], remove_indices, axis=0)

        # go through each graph name
        for graph_name in self.graph_names:
            # add/remove vertices from the graph
            self.__dict__[graph_name].add_vertices(num_added)
            self.__dict__[graph_name].delete_vertices(remove_indices)

        # change total number of agents and print info to terminal
        self.number_agents += num_added
        # print("\tAdded " + str(num_added) + " agents")
        # print("\tRemoved " + str(num_removed) + " agents")

        # clear the hatching/removing arrays for the next step
        self.hatching[:] = False
        self.removing[:] = False

    def cell_fate(self):

        # Bernoulli Probabilities:

        # cho_types = np.argwhere(self.cell_type == 0)
        # if len(cho_types) > 0:
        #     committed_cells = np.random.binomial(1, self.transition_rate, len(cho_types))
        #     aba_cells = np.random.binomial(1, self.aba_ratio / (self.aba_ratio + self.dox_ratio), sum(committed_cells)) + 1
        #     cell_fate = aba_cells
        #     start = np.where(self.cell_type == 0)[0][0]
        #     self.cell_type[start:start+len(cell_fate)] = cell_fate
        # self.update_colors()

        transition_threshold = calculate_transition_rate(self.transition_rate, self.current_step)
        committed_cells = self.transition < transition_threshold
        for i in range(self.number_agents):
            if self.cell_type[i] == 0 and committed_cells[i]:
                self.cell_type[i] = np.random.binomial(1, self.aba_ratio / (self.dox_ratio + self.aba_ratio),  1) + 1
        self.update_colors()

    def update_colors(self):
        ref = np.array([self.cho_color, self.dox_color, self.aba_color])
        self.colors = ref[self.cell_type]

    @record_time
    def move_parallel(self):
        edges = np.asarray(self.neighbor_graph.get_edgelist())
        num_edges = len(edges)
        edge_forces = np.zeros((num_edges, 2, 3))
        neighbor_forces = np.zeros((self.number_agents, 3))
        # get adhesive/repulsive forces from neighbors and gravity forces
        edge_forces = get_neighbor_forces(num_edges, edges, edge_forces, self.locations, self.center, self.cell_type,
                                          self.cell_rad, u_bb=self.u_bb, u_rb=self.u_rb, u_rr=self.u_rr, u_yb=self.u_yb,
                                          u_ry=self.u_ry, u_yy=self.u_yy, alpha=0, u_repulsion=self.u_repulsion)
        neighbor_forces = convert_edge_forces(num_edges, edges, edge_forces, neighbor_forces)
        noise_vector = np.ones((self.number_agents, 3)) * self.alpha * (2 * np.random.rand(self.number_agents, 3) - 1)
        neighbor_forces = neighbor_forces + noise_vector
        if self.gravity > 0:
            net_forces = np.zeros((self.number_agents, 3))
            gravity_forces = get_gravity_forces(self.number_agents, self.locations, self.center,
                                                325, net_forces, grav=self.gravity)
            neighbor_forces = neighbor_forces + gravity_forces
        for i in range(self.number_agents):
            vec = neighbor_forces[i]
            sum = vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2
            if sum != 0:
                neighbor_forces[i] = neighbor_forces[i] / (sum ** (1/2))
            else:
                neighbor_forces[i] = 0
        # update locations based on forces
        self.locations += 2 * self.velocity * self.cell_rad * neighbor_forces
        # check that the new location is within the space, otherwise use boundary values
        self.locations = np.where(self.locations > self.well_rad, self.well_rad, self.locations)
        self.locations = np.where(self.locations < 0, 0, self.locations)

    #Unused..
    @record_time
    def reproduce(self, ts):
        """ If the agent meets criteria, hatch a new agent.
        """
        # increase division counter by time step for all agents
        self.division_set += ts

        #go through all agents marking for division if over the threshold
        if self.replication_type == 'Contact_Inhibition':
            adjacency_matrix = self.neighbor_graph.get_adjacency()
            for index in range(self.number_agents):
                if self.division_set[index] > self.div_thresh[index]:
                    # 12 is the maximum number of cells that can surround a cell
                    if np.sum(adjacency_matrix[index,:]) < 12:
                        self.mark_to_hatch(index)
        if self.replication_type == 'Default':
            for index in range(self.number_agents):
                if self.division_set[index] > self.div_thresh[index]:
                    self.mark_to_hatch(index)
        if self.replication_type == 'None':
            return

    @classmethod
    def simulation_mode_0(cls, name, output_dir, yaml_file="general.yaml"):
        """ Creates a new brand new simulation and runs it through
            all defined steps.
        """
        # make simulation instance, update name, and add paths
        sim = cls(yaml_file)
        sim.name = name
        sim.set_paths(output_dir)

        # set up the simulation agents and run the simulation
        sim.full_setup()
        sim.run_simulation()

    def save_params(self, path):
        """ Add the instance variables to the Simulation object based
            on the keys and values from a YAML file.
        """
        # load the dictionary
        params = template_params(path)

        # iterate through the keys adding each instance variable
        with open(self.main_path + "parameters.txt", "w") as parameters:
            for key in list(params.keys()):
                parameters.write(f"{key}: {params[key]}\n")
        parameters.close()

    @classmethod
    def start_sweep(cls, output_dir, yaml_file, name, mode):
        """ Configures/runs the model based on the specified
            simulation mode.
        """
        # check that the output directory exists and get the name/mode for the simulation
        output_dir = backend.check_output_dir(output_dir)

        # new simulation
        if mode == 0:
            # first check that new simulation can be made and run that mode.
            # i = 0
            # while os.path.isdir(output_dir + name):
            #     name = name + f'_{i}'
            #     i += 1
            print(f'Starting {name} ...')
            name = backend.check_existing(name, output_dir, new_simulation=True)
            cls.simulation_mode_0(name, output_dir, yaml_file=yaml_file)

        # existing simulation
        else:
            # check that previous simulation exists
            name = backend.check_existing(name, output_dir, new_simulation=False)

            # call the corresponding mode
            if mode == 1:
                cls.simulation_mode_1(name, output_dir)    # continuation
            elif mode == 2:
                cls.simulation_mode_2(name, output_dir)    # images to video
            elif mode == 3:
                cls.simulation_mode_3(name, output_dir)    # archive simulation
            else:
                raise Exception("Mode does not exist!")


#EDIT WHICH YAML FILE USED IN simulation_mode_0.
if __name__ == "__main__":
    TestSimulation.start("/Users/andrew/PycharmProjects/CHO_adhesion_model/outputs/")