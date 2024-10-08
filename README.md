# Web based 2D Wind Tunnel Simulator
#### Video Demo:  <URL HERE>
#### Description: The Web-based 2D Wind Tunnel Simulator is an interactive educational tool designed to simulate airflow around a simplified car (half circle) in a two-dimensional environment. This project aims to provide students, engineers, and aerodynamics enthusiasts with a user-friendly platform to visualize and understand fluid dynamics principles without the need for complex physical setups or expensive software.

## Key features of the simulator include:
1. Visualization of airflow patterns
2. Adjustable wind speed and car size
3. Interactive flow field display

</div>

## Structure
The 2D Wind Tunnel Simulator includes the following files:
* `README`: this file
* `LICENSE`: the MIT License file
* `index.html`: main page for Simulator User Interface with control buttons, visualization options and simulation parameters
* `styles.css`: the file for Styling and User Interface Design 
* `windtunnel.js`: the javascript file for Fluid Solver, Visualization and Interactive Control
  
</div>

## Javascript Implementation in `windtunnel.js`
* `Class Fluid`: The heart of 2D Wind Tunnel Simulator, which implements the core fluid dynamics simulation. 
  * `constructor(numX, numY, h, density)`: Initializes the fluid grid (with extra boundary cells and properties (velocity fields, pressure, smoke density).
  * `solveIncomp(dt, maxIters)`: Solves for incompressibility using Jacobi iteration. Calculates pressure and applies pressure corrections to enforce the divergence-free condition.
  * `fieldCalc(x, y, field)`:
  * `fieldCalc(x, y, field)`:
  * `fieldCalc(x, y, field)`: 
## 
Performs bilinear interpolation to calculate field values (velocity or smoke density) at arbitrary positions.
## avgU(i, j) and avgV(i, j)
Calculate average velocities at cell centers, crucial for advection steps.
## applyBoundaryConditions()
Handles boundary conditions by extrapolating velocities to boundary cells.
## advectVel(dt)
Implements semi-Lagrangian advection for the velocity field.
## advectSmoke(dt)
Advects the smoke density field using the semi-Lagrangian method.
## simulate(dt, maxIters)
Main simulation step that combines all components: pressure solving, boundary condition application, and advection.
# Design Choices and Rationale
## Semi-Lagrangian Advection
We chose the semi-Lagrangian method for advection due to its unconditional stability, allowing for larger time steps without numerical instability.
This method provides a good balance between accuracy and computational efficiency, crucial for real-time web-based simulations.
## Jacobi Iteration for Pressure Solving
The Jacobi method is used for its simplicity and ease of implementation.
While not the fastest method, it's sufficiently efficient for our 2D simulation and easier to understand, making the code more accessible for educational purposes.
## Staggered Grid
We use a staggered grid arrangement where pressure is stored at cell centers, and velocity components are stored on cell faces.
This choice helps in reducing numerical instabilities and provides more accurate pressure gradient calculations.
## Over-Relaxation Option
The inclusion of an over-relaxation parameter (controlled by scene.overRelaxation) allows for faster convergence in the pressure solving step, at the cost of potential instability if set too high.
## Smoke Density Simulation
Adding smoke density simulation enhances the visual feedback of the fluid flow, making the simulation more intuitive and engaging for users.
## Bilinear Interpolation
The fieldCalc method uses bilinear interpolation for smooth and accurate sampling of field values between grid points.
## Boundary Handling
Special care is taken to handle boundaries correctly, ensuring the simulation remains stable and realistic at the edges of the domain.
These design choices aim to create a fluid simulation that is computationally efficient enough to run in real-time in a web browser while still providing accurate and visually appealing results. The code structure also allows for easy extension and modification, making it suitable for educational purposes and further experimentation.

# Scene Setup and Initialization
The setupScene() function plays a vital role in setting up the initial conditions for our 2D Wind Tunnel Simulator. This function initializes the fluid simulation domain, sets boundary conditions, and configures various simulation parameters.
Key Components and Functionality
## Resolution Setting
Allows switching between high (240x240) and standard (80x80) resolution.
Adjusts the time step for high-resolution mode to maintain stability.
## Domain Configuration
Sets up the simulation domain based on the chosen resolution.
Calculates the number of cells in X and Y directions to maintain the aspect ratio.
## Fluid Initialization
Creates a new Fluid instance with the calculated parameters.
## Boundary Conditions
Sets up solid boundaries (s = 0) at the top, bottom, and left edges of the domain.
Initializes the inlet velocity at the left boundary.
## Smoke Initialization
Initializes a vertical strip of smoke near the left boundary for flow visualization.
## Obstacle Placement
Places the obstacle (simulated car) in the wind tunnel.
## UI Synchronization
Ensures that the UI controls reflect the current simulation settings.
# Design Choices and Rationale
## Flexible Resolution:
The option to switch between high and standard resolution allows users to balance between simulation accuracy and performance based on their hardware capabilities.
## Aspect Ratio Preservation:
The domain width is calculated to maintain the correct aspect ratio, ensuring the simulation accurately represents the intended physical space.
## Density Setting:
Using a realistic fluid density (998.0, close to water) provides physically accurate behavior, though the simulation is scale-invariant.
## Boundary Conditions:
Solid boundaries are set up to create a confined wind tunnel effect.
The inlet velocity is set at the left boundary to simulate wind.
## Smoke Initialization:
A thin strip of smoke is initialized for better flow visualization, enhancing the educational value of the simulation.
## UI Synchronization:
Ensures that the user interface accurately reflects the current state of the simulation, improving user experience and preventing confusion.
This setup function demonstrates the careful balance between physical accuracy, computational efficiency, and user experience in creating an educational fluid dynamics simulation. It provides a solid foundation for the simulation while allowing for easy adjustments and extensions.


# Color Management and Visualization
Our 2D Wind Tunnel Simulator employs sophisticated color management techniques to provide clear and intuitive visual representations of the fluid dynamics data. Two key functions handle this aspect of the simulation:
## setColor(r, g, b)
Converts normalized RGB values (0-1) to the 0-255 range.
Sets both the fill and stroke styles of the canvas context to the specified color.
Design Choice:
Using a single function to set both fill and stroke styles ensures consistency in color application across different drawing operations.
The function accepts normalized RGB values, making it easier to work with color calculations in the simulation logic.
## getSciColor(val, minVal, maxVal)
Implements a scientific color map for visualizing scalar fields (like pressure or velocity magnitude).
Maps a value within a given range to a color using a four-segment color gradient.
## Color Mapping Process:
Normalizes the input value to a range of [0,1].
Divides the color spectrum into four segments.
Determines which segment the normalized value falls into.
Interpolates the color within that segment.
## Color Scheme:
Blue (0.0) → Cyan (0.25) → Green (0.5) → Yellow (0.75) → Red (1.0)
Design Choices:
Four-Segment Gradient:
Provides a wide range of distinct colors, making it easier to discern small variations in the data.
The chosen color progression is intuitive for most users (cool colors for low values, warm colors for high values).
Normalization:
Allows the function to work with any range of input values, making it versatile for different types of data (e.g., pressure, velocity).
Clamping:
Values outside the specified range are clamped to the nearest endpoint, preventing unexpected color outputs.
RGBA Output:
Returns color as an array of RGBA values (0-255 range), making it compatible with various drawing functions and allowing for alpha channel manipulation if needed.
# Importance in the Simulation
These color functions play a crucial role in making the simulation data visually comprehensible:
Intuitive Data Representation: The scientific color map allows users to quickly understand the distribution and intensity of various fluid properties across the simulation domain.
Enhanced Visual Feedback: By mapping data to colors, subtle changes in the fluid behavior become more apparent, enhancing the educational value of the simulation.
Flexibility: The functions are designed to work with different types of data and visualization methods (e.g., pressure fields, velocity magnitudes, streamlines), providing a consistent visual language throughout the simulation.
Performance Consideration: By pre-calculating colors based on data values, we can efficiently render large amounts of visual information without recalculating colors for each frame or pixel.
These color management functions contribute significantly to the user experience and educational effectiveness of the 2D Wind Tunnel Simulator, turning complex numerical data into intuitive visual representations.

# Coordinate Conversion
In our 2D Wind Tunnel Simulator, we need to map the simulation's coordinate system to the canvas coordinate system for proper rendering. Two simple yet essential functions handle this conversion:
cX(x)
Converts the x-coordinate from simulation space to canvas space.
Scales the x-coordinate by cScale (canvas scale factor).
cY(y)
Converts the y-coordinate from simulation space to canvas space.
Scales the y-coordinate and inverts it (subtracts from canvas height).
# Importance in the Simulation
Accurate Visualization: Ensures that the fluid dynamics simulation is accurately represented on the canvas, maintaining proper proportions and orientations.
Flexibility: Allows for easy adjustments to the view (like zooming or resizing) by modifying the cScale factor without changing the underlying simulation code.
Abstraction: Provides a clean interface between the simulation's coordinate system and the canvas rendering, making it easier to modify either part independently.
Consistency: Ensures that all elements of the simulation (fluid particles, obstacles, boundaries) are consistently positioned relative to each other when rendered.
These coordinate conversion functions, while simple, play a crucial role in bridging the gap between the mathematical model of the fluid simulation and its visual representation on the screen. They contribute significantly to the accuracy and intuitiveness of the visual output, enhancing the educational value of the 2D Wind Tunnel Simulator.

# Visualization and Rendering
The map() function is the core rendering function of our 2D Wind Tunnel Simulator. It's responsible for visualizing the fluid simulation data on the canvas, including pressure, smoke, and streamlines.
# Key Components and Functionality
## Canvas Preparation
Clears the canvas before each render to prevent ghosting.
## Resolution Adjustment
Adjusts cell size based on the chosen resolution for optimal visualization.
## Pressure Range Calculation
Determines the pressure range for color mapping.
## Pixel-by-Pixel Rendering
Uses getImageData and putImageData for efficient pixel manipulation.
Allows for fine-grained control over the visualization.
## Visualization Modes
Supports multiple visualization modes:
Pressure visualization with color mapping
Smoke density visualization
Combined pressure and smoke visualization
## Streamline Rendering
Draws streamlines to visualize fluid flow patterns.
Adjusts streamline density based on resolution.
## Pressure Range Display
Displays the current pressure range when pressure visualization is active.

# Obstacle Handling
In our 2D Wind Tunnel Simulator, the obstacle represents a simplified car model (a half-circle). Two key functions manage the obstacle's rendering and physics interaction:
## drawObstacle()
This function is responsible for rendering the obstacle on the canvas.
## Key Features:
Draws a half-circle shape to represent the car.
Adjusts the appearance based on the current visualization mode (pressure or normal).
Adds a bottom line to complete the car silhouette.

## setObstacle(x, y, reset)
This function handles the positioning and physics interaction of the obstacle within the fluid simulation.
## Key Features:
Updates the obstacle's position in the simulation.
Calculates the obstacle's velocity if it's being moved.
Modifies the fluid simulation grid to account for the obstacle's presence.

These functions play a crucial role in making the 2D Wind Tunnel Simulator an interactive and educational tool. By allowing users to visualize and interact with a simplified car model in a fluid environment, it provides intuitive insights into basic aerodynamics principles. The careful design of these functions balances visual clarity, physical accuracy, and computational efficiency, enhancing the overall user experience and educational value of the simulator.

# User Interaction and Simulation Control
Our 2D Wind Tunnel Simulator incorporates several functions to handle user interactions and control the simulation flow. These functions enhance the interactivity and usability of the simulator, making it a powerful educational tool.
## Mouse Interaction
function startDrag(x, y) { /* ... */ }
function drag(x, y) { /* ... */ }
function endDrag() { /* ... */ }
## Simulation Control
function togglePause() { /* ... */ }
function stepMotion() { /* ... */ }
## Main Simulation Loop
function simulate() { /* ... */ }
function update() { /* ... */ }

# Initialization and UI Setup
The init() function plays a vital role in setting up the 2D Wind Tunnel Simulator. It initializes the simulation scene and configures the user interface controls, ensuring that the simulator is ready for user interaction.
# Design Choices and Rationale
## Modular Initialization
Separating the initialization into setupScene() and init() allows for easier maintenance and potential reuse of setup logic.
## Real-time UI Responsiveness
Event listeners are set up to respond immediately to user input, providing instant feedback and a more engaging user experience.
## Direct Simulation Updates
Changes in UI controls (like wind velocity) directly modify the simulation parameters, ensuring tight integration between the UI and the underlying physics model.
## Resolution Switching
The ability to switch between standard and high-resolution modes allows users to balance between performance and detail based on their hardware capabilities.
## Error Handling
The use of a try-catch block ensures that any errors during initialization are caught and logged, improving the robustness of the application.

The init() function is fundamental to creating an interactive and user-friendly experience in the 2D Wind Tunnel Simulator. It sets the stage for user exploration and learning by ensuring that all components of the simulator – from the underlying physics model to the user interface controls – are properly initialized and interconnected. This careful setup is key to making complex fluid dynamics concepts accessible and engaging to users, regardless of their technical background.


TODO List:
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
11. Implement Hi/Low resolution switch
12. Performance Optimization
13. Unit Test
14. Publish project to GitHub Pages or similar hosting service
15. Create project showcase video and update README with link

Future Development:
1. Complicated Car-shape obstacle (obstacle detector with arbitary shapes)
2. User-defined obstacle

References:
1. [Chentanez, N., Muller, M. (2011). NVIDIA PhysX Research, Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid](https://matthias-research.github.io/pages/publications/tallCells.pdf#:~:text=We%20present%20a%20new%20Eulerian%20fluid)
2. [Chentanez, N., Muller, M. (2011). NVIDIA, Real-Time Eulerian Water Simulation Using a Restricted Tall Cell Grid, SIGGRAPH2011](https://matthias-research.github.io/pages/publications/tallCellsSlides.pdf#:~:text=%E2%80%A2%20GPU%20friendly%20tall%20cell%20grid)
3. [Ash, M., (2006). Fluid Simulation for Dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html#:~:text=With%20that%20goal%20in%20mind,%20I'm#:~:text=With%20that%20goal%20in%20mind,%20I'm)
4. [NoBS Code, (2024). The Midpoint Circle Algorithm Explained Step by Step](https://www.youtube.com/watch?v=hpiILbMkF9w)

