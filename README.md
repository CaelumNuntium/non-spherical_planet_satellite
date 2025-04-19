# Modeling the motion of a satellite of a non-spherical planet

The program numerically simulates the motion of a close satellite of a non-spherical planet with given value of $J_2$ harmonic in the expansion of the gravitational potential. The program allows you to compare the result of numerical modeling with the analytical model of the motion of a close satellite of a non-spherical planet, obtained by averaging the Laplace equations.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/CaelumNuntium/non-spherical_planet_satellite.git
cd non-spherical_planet_satellite
```

2. Build program:
```bash
make
```
## Usage

1. First, set the parameters in *config.ini*

2. Run the program:

    * On Linux:
    ```bash
    ./nsps
    ```
    
    * On Windows:
    ```cmd
    nsps.exe
    ```

3. Visualization:
```bash
python plot.py
```
Files *semimajor_axis.png*, *eccentricity.png*, *inclination.png*, *asc_node.png*, *pericenter_argument.png*, *mean_anomaly.png* will contain graphs of the osculating elements $a$, $e$, $i$, $\Omega$, $g$, $M$ of the satellite's orbit in the numerical and analytical models of motion. Files *delta_a.png*, *delta_e.png*, *delta_i.png*, *delta_Omega.png*, *delta_g.png*, *delta_M.png* will contain graphs of the differences between the models.
