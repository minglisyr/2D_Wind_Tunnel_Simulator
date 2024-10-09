# Web based 2D Wind Tunnel Simulator
#### Video Demo:  <URL HERE>
[Project Web Hosting](https://minglisyr.github.io/2D_Wind_Tunnel_Simulator)
#### Description: The Web-based 2D Wind Tunnel Simulator is an interactive educational tool designed to simulate airflow around a simplified car (half circle) in a two-dimensional environment. This project aims to provide students, engineers, and aerodynamics enthusiasts with a user-friendly platform to visualize and understand fluid dynamics principles without the need for complex physical setups or expensive software.

![Demo](https://github.com/minglisyr/2D_Wind_Tunnel_Simulator/blob/main/Demo.gif)

## Key features of the simulator include:
1. Visualization of airflow patterns
2. Adjustable wind speed and car size
3. Interactive flow field display

</div>

## Structure
The 2D Wind Tunnel Simulator includes the following files:
* `README.md`: this file
* `LICENSE`: the MIT License file
* `index.html`: the Simulator User Interface with control buttons, visualization options and simulation parameters
* `styles.css`: the file for Styling and User Interface Design 
* `windtunnel.js`: the javascript file for Fluid Solver, Visualization and Interactive Control
  
</div>

## Javascript Implementation in `windtunnel.js`
  #### Core Algorithm of Fluid Simulator
  * `Class Fluid`: The heart of 2D Wind Tunnel Simulator, which implements the core fluid dynamics simulation. 
    * `constructor(numX, numY, h)`: Initializes the fluid grid (with extra boundary cells) and properties (velocity fields, pressure, smoke density).
    * `solveIncomp(dt, maxIters)`: Solves for incompressibility using Jacobi iteration method. Calculates pressure and applies pressure corrections to enforce the divergence-free condition.
    * `fieldCalc(x, y, field)`: Performs bilinear interpolation to calculate field values (velocity or smoke density) at arbitrary positions.
    * `applyBoundaryConditions()`: Handles boundary conditions by extrapolating velocities to boundary cells.
    * `advectVel(dt) and advectSmoke(dt)`: Advects the velocity field and smoke density field using the semi-Lagrangian method.
    * `main(dt, maxIters)`: Main simulation step that combines all components: pressure solving, boundary condition application, and advection.

    > [!note]
    > * The `Semi-Lagrangian Method` was chosen for advection due to its unconditional stability, allowing for larger time steps without numerical instability.
    > * The `Jacobi Iteration Method` was chosen for pressure solving for its simplicity and ease of implementation. While not the fastest method, it's sufficiently efficient for our 2D simulation and easier to understand, making the code more accessible for educational purposes. For higher resolution real-time solving, `Multigrid solver` might be a better choice.
    > * The `Staggered Grid` is used to store pressure at cell centers, and velocity components are stored on cell faces, which helps in reducing numerical instabilities and provides more accurate pressure gradient calculations.
    > * The `Over-Relaxation` is used for faster convergence in the pressure solving step, at the cost of potential instability if set too high.

  #### Setup Scence and Initialization
  * `setupScene()`: Initializes the fluid simulation domain, sets boundary conditions, and configures various simulation parameters.
    * `Resolution Setting`: Allows switching between high (240x240) and standard (80x80) resolution. Adjusts the time step for high-resolution mode to maintain stability.
    * `Domain Configuration`: Sets up the simulation domain based on the chosen resolution. Calculates the number of cells in X and Y directions to maintain the aspect ratio.
    * `Fluid Initialization`: Creates a new Fluid instance with the calculated parameters.
    * `Boundary Conditions`: Sets up solid boundaries (s = 0) at the top, bottom, and left edges of the domain. Initializes the inlet velocity at the left boundary.
    * `Smoke Initialization`: Initializes a vertical strip of smoke near the left boundary for flow visualization.
    * `Obstacle Placement`: Places the obstacle (simulated car) in the wind tunnel.
    * `UI Synchronization`: Ensures that the UI controls reflect the current simulation settings.

  #### Color Management and Visualization
  
  * `getSciColor(val, minVal, maxVal)`: Implements a scientific color map for visualizing scalar fields (like pressure or velocity magnitude).
    * `Color Mapping Process`:  Normalize the input value to a range of [0,1] → Divide the color spectrum into four segments → Determine which segment the normalized value falls into → Interpolate the color within that segment
    * `Color Scheme`: $${\color{blue}Blue}$$ (0.0) → $${\color{cyan}Cyan}$$ (0.25) → $${\color{green}Green}$$ (0.5) → $${\color{yellow}Yellow}$$ (0.75) → $${\color{red}Red}$$ (1.0)

  #### Coordinate Conversion
  * `cX(x)`: Converts the x-coordinate from simulation space to canvas space. Scales the x-coordinate by cScale (canvas scale factor).
  * `cY(y)`: Converts the y-coordinate from simulation space to canvas space. Scales the y-coordinate and inverts it (subtracts from canvas height).

  #### Visualization and Rendering
  * `map()`: Core rendering function, which is responsible for visualizing the fluid simulation data on the canvas, including pressure, smoke, and streamlines.
    * `Canvas Preparation`: Clears the canvas before each render to prevent ghosting.
    * `Resolution Adjustment`: Adjusts cell size based on the chosen resolution for optimal visualization.
    * `Pressure Range Calculation`: Determines the pressure range for color mapping.
    * `Pixel-by-Pixel Rendering`: Uses getImageData and putImageData for efficient pixel manipulation. Allows for fine-grained control over the visualization.
    * `Visualization Modes`: Supports multiple visualization modes: 1) Pressure 2) Smoke Density 3) Combined
    * `Streamline Rendering`: Draws streamlines to visualize fluid flow patterns. Adjusts streamline density based on resolution.
    * `Pressure Range Display`: Displays the current pressure range when pressure visualization is active.

  #### Obstacle Handling
  * `setObstacle(x, y, reset)`: This function handles the positioning and physics interaction of the obstacle within the fluid simulation.
    * Updates the obstacle's position in the simulation.
    * Calculates the obstacle's velocity if it's being moved.
    * Modifies the fluid simulation grid to account for the obstacle's presence.
  * `drawObstacle()`: This function is responsible for rendering the obstacle on the canvas.
    * Draws a half-circle shape to represent the car.
    * Adjusts the appearance based on the current visualization mode (pressure or normal).
    * Adds a bottom line to complete the car silhouette.

  #### User Interaction and Simulation Control
| Mouse Interaction         | Simulation Control     | Main Simulation Loop |
|--------------|-----------|------------|
| startDrag(x, y) | togglePause( )      | simulate( )        |
| drag(x, y)      | stepMotion( )  | update( )       |
| endDrag( )       |   |        |

  #### Initialization and UI Setup
  * `init()`: Initializes the simulation scene and configures the user interface controls, ensuring that the simulator is ready for user interaction.
    * `Modular Initialization`: Separating the initialization into setupScene() and init() allows for easier maintenance and potential reuse of setup logic.
    * `Real-time UI Responsiveness`: Event listeners are set up to respond immediately to user input, providing instant feedback and a more engaging user experience.
    * `Direct Simulation Updates`: Changes in UI controls (like wind velocity) directly modify the simulation parameters, ensuring tight integration between the UI and the underlying physics model.
    * `Resolution Switching`: The ability to switch between standard and high-resolution modes allows users to balance between performance and detail based on their hardware capabilities.
    * `Error Handling`: The use of a try-catch block ensures that any errors during initialization are caught and logged, improving the robustness of the application.
    ##### [!NOTE]
    The init() function sets the stage for user exploration and learning by ensuring that all components of the simulator – from the underlying physics model to the user interface controls – are properly initialized and interconnected. This careful setup is key to making complex fluid dynamics concepts accessible and engaging to users, regardless of their technical background.


</div>

## TODO List (for CS50x):
1. Create basic file structure (HTML, CSS, JavaScript)
2. Design main layout
3. Create responsive grid for simulation area
4. Core Simulation Logic
5. Implement color-coded velocity and pressure mapping
6. Implement streamline visualization option
7. Implement smoke visualization option
8. User Controls (Restart, Pause/Resume, Step Forward)
9. User interactive with obstacle positioning and movement
10. Implement wind speed control and car size control
11. Implement High/Standard resolution switch
12. Performance Optimization
13. Unit Test
14. Publish project to GitHub Pages or similar hosting service
15. Create project showcase video and update README with link

</div>

## Future Development:
1. Complicated/User-defined car-shape obstacle (obstacle detector with arbitary shapes)
2. Multigrid solver for higher resolution implementation

</div>

## References:
1. [Chentanez, N., Muller, M. (2011). NVIDIA PhysX Research, Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid](https://matthias-research.github.io/pages/publications/tallCells.pdf#:~:text=We%20present%20a%20new%20Eulerian%20fluid)
2. [Chentanez, N., Muller, M. (2011). NVIDIA, Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid, SIGGRAPH2011](https://matthias-research.github.io/pages/publications/tallCellsSlides.pdf#:~:text=%E2%80%A2%20GPU%20friendly%20tall%20cell%20grid)
3. [Ash, M., (2006). Fluid Simulation for Dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html#:~:text=With%20that%20goal%20in%20mind,%20I'm#:~:text=With%20that%20goal%20in%20mind,%20I'm)
4. [NoBS Code, (2024). The Midpoint Circle Algorithm Explained Step by Step](https://www.youtube.com/watch?v=hpiILbMkF9w)

