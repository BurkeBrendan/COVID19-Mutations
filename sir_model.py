import random
from cellular_automaton import *


RECOVERED = [1.0]
INFECTED = [0.1]
SUSCEPTIBLE = [0]

IMMUNE = [3.0]

MUTANT_INFECTED = [1.1]
MUTANT_RECOVERED = [2.0]

init_mode = 0


class SIRModel(Rule):
    random_seed = random.seed(13)
    cycle = 0

    def init_state(self, cell_coordinate):
        if init_mode == 0:
            if cell_coordinate == (20, 20):
                init = 0.1
            else:
                init = 0.0
        else:
            rand = random.randrange(0, 80, 1)
            init = max(.0, float(rand - 78))
            init = init / 10.0
        return [init]

    @staticmethod
    def __infect():
        return random.random() < 0.25

    def evolve_cell(self, last_cell_state, neighbors_last_states):
        self.cycle = self.cycle + 1
        if self.cycle == 600000:
            return MUTANT_INFECTED
        new_cell_state = last_cell_state
        sick_neighbours = self.__count_sick_neighbours(neighbors_last_states)
        mutated_sick_neighbors = self.__count_mutated_sick_neighbours(neighbors_last_states)
        prob_infected = self.__infect()
        if last_cell_state == IMMUNE:
            new_cell_state = IMMUNE
        elif (last_cell_state == RECOVERED or last_cell_state == SUSCEPTIBLE or last_cell_state == INFECTED) and mutated_sick_neighbors > 0 and prob_infected:
            new_cell_state = MUTANT_INFECTED
        elif last_cell_state == RECOVERED:
            new_cell_state == RECOVERED
        elif last_cell_state[0] > 1.0 and last_cell_state[0] < 1.9:
            new_cell_state = [last_cell_state[0] + 0.1]
        elif last_cell_state[0] >= 1.9 and last_cell_state[0] <= 2.0:
            new_cell_state = MUTANT_RECOVERED
        elif last_cell_state == SUSCEPTIBLE and sick_neighbours > 0 and prob_infected:
            new_cell_state = INFECTED
        elif last_cell_state == SUSCEPTIBLE and (sick_neighbours == 0 or not prob_infected):
            new_cell_state = SUSCEPTIBLE
        elif last_cell_state[0] < 0.9:
            new_cell_state = [last_cell_state[0] + 0.1]
        elif last_cell_state[0] >= 0.9 and last_cell_state[0] <= 1.0:
            if random.random() < 0.4:
                new_cell_state = IMMUNE
            else:
                new_cell_state = RECOVERED
        else:
            new_cell_state = MUTANT_RECOVERED
        return new_cell_state

    @staticmethod
    def __count_sick_neighbours(neighbours):
        an = []
        for n in neighbours:
            if n[0] > 0.0 and n[0] < 1.0:
                an.append(1)
        return len(an)

    @staticmethod
    def __count_mutated_sick_neighbours(neighbours):
        an = []
        for n in neighbours:
            if n[0] > 1.0 and n[0] < 2.0:
                an.append(1)
        return len(an)

    def get_state_draw_color(self, current_state):
        color = [100, 100, 255]
        if current_state[0] == 0.0:
            color = [0, 255, 0]
        elif current_state[0] < 1.0:
            color = [255, 0, 0]
        elif current_state[0] > 1.0 and current_state[0] < 2.0:
            color = [255, 165, 0]
        elif current_state[0] == 2.0:
            color = [0, 0, 0]
        elif current_state[0] == 3.0:
            color = [100, 100, 100]
        return color


if __name__ == "__main__":
    neighborhood = MooreNeighborhood(EdgeRule.IGNORE_MISSING_NEIGHBORS_OF_EDGE_CELLS)
    ca = CAFactory.make_multi_process_cellular_automaton(dimension=[250, 250],
                                                         neighborhood=neighborhood,
                                                         rule=SIRModel,
                                                         processes=1)
    ca_window = CAWindow(cellular_automaton=ca, evolution_steps_per_draw=1)
