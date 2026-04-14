from simulate import BalloonSimulation
from plotting import run_all_plots


def main():
    sim = BalloonSimulation()
    sim.run()
    run_all_plots(sim.sim_state, sim)


if __name__ == "__main__":
    main()