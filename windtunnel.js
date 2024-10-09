const canvas = document.getElementById("myCanvas");
const ctx = canvas.getContext("2d", {
    willReadFrequently: true
});

canvas.width = window.innerWidth - 20;
canvas.height = window.innerHeight - 100;
canvas.focus();

const simHeight = 1.005;
const cScale = canvas.height / simHeight;
const simWidth = canvas.width / cScale;

const U_FIELD = 1;
const V_FIELD = 2;
const S_FIELD = 9;

// ----------------- Fluid Core Algorithm ------------------------------
class Fluid {
    constructor(numX, numY, h) {
        // Initialize fluid properties
        this.numX = numX + 2; // Grid cells in X direction (including 2 boundary cells)
        this.numY = numY + 2; // Grid cells in Y direction (including 2 boundary cells)
        this.numCells = this.numX * this.numY; // Total number of cells
        this.h = h; // Cell size

        // Initialize arrays for fluid properties
        this.u = new Float32Array(this.numCells).fill(0.0);
        this.v = new Float32Array(this.numCells).fill(0.0);
        this.newU = new Float32Array(this.numCells).fill(0.0);
        this.newV = new Float32Array(this.numCells).fill(0.0);

        this.p = new Float32Array(this.numCells).fill(0.0);
        this.s = new Float32Array(this.numCells).fill(0.0);
        this.smo = new Float32Array(this.numCells).fill(1.0); // Initialize smoke density to 1.0 everywhere
        this.newSMO = new Float32Array(this.numCells).fill(1.0);
    }

    solveIncomp(dt, maxIters) {
        const n = this.numY;
        const cp = this.h / dt; // Pressure coefficient

        // Jacobian iteration to solve for pressure and enforce incompressibility
        for (let iter = 0; iter < maxIters; iter++) {
            for (let i = 1; i < this.numX - 1; i++) {
                for (let j = 1; j < this.numY - 1; j++) {
                    if (this.s[i * n + j] == 0.0) continue; // Skip non-fluid cells

                    // Calculate divergence
                    const sx0 = this.s[(i - 1) * n + j];
                    const sx1 = this.s[(i + 1) * n + j];
                    const sy0 = this.s[i * n + j - 1];
                    const sy1 = this.s[i * n + j + 1];
                    const s = sx0 + sx1 + sy0 + sy1;
                    if (s == 0.0) continue;
                    const div = this.u[(i + 1) * n + j] - this.u[i * n + j] + this.v[i * n + j + 1] - this.v[i * n + j];

                    // Calculate pressure correction
                    let p = -div / s;
                    p *= scene.overRelaxation;
                    this.p[i * n + j] += cp * p;

                    // Apply pressure correction to velocities
                    this.u[i * n + j] -= sx0 * p;
                    this.u[(i + 1) * n + j] += sx1 * p;
                    this.v[i * n + j] -= sy0 * p;
                    this.v[i * n + j + 1] += sy1 * p;
                }
            }
        }
    }

    fieldCalc(x, y, field) {
        // Bilinear interpolation to field calculation at given (x, y) position
        const n = this.numY;
        const h = this.h;

        x = Math.max(Math.min(x, this.numX * h), h);
        y = Math.max(Math.min(y, this.numY * h), h);

        let dx = 0.0;
        let dy = 0.0;
        let f;

        switch (field) {
            case U_FIELD:
                f = this.u;
                dy = 0.5 * h;
                break;
            case V_FIELD:
                f = this.v;
                dx = 0.5 * h;
                break;
            case S_FIELD:
                f = this.smo;
                dx = 0.5 * h;
                dy = 0.5 * h;
                break;
        }

        const x0 = Math.min(Math.floor((x - dx) / h), this.numX - 1);
        const tx = ((x - dx) - x0 * h) / h;
        const x1 = Math.min(x0 + 1, this.numX - 1);
        const y0 = Math.min(Math.floor((y - dy) / h), this.numY - 1);
        const ty = ((y - dy) - y0 * h) / h;
        const y1 = Math.min(y0 + 1, this.numY - 1);
        const sx = 1.0 - tx;
        const sy = 1.0 - ty;

        return sx * sy * f[x0 * n + y0] +
            tx * sy * f[x1 * n + y0] +
            tx * ty * f[x1 * n + y1] +
            sx * ty * f[x0 * n + y1];
    }

    applyBoundaryConditions() {
        const n = this.numY;

        // Extrapolate horizontal velocities to boundary cells
        for (let i = 0; i < this.numX; i++) {
            this.u[i * n + 0] = this.u[i * n + 1];
            this.u[i * n + this.numY - 1] = this.u[i * n + this.numY - 2];
        }

        // Extrapolate vertical velocities to boundary cells
        for (let j = 0; j < this.numY; j++) {
            this.v[0 * n + j] = this.v[1 * n + j];
            this.v[(this.numX - 1) * n + j] = this.v[(this.numX - 2) * n + j];
        }
    }

    advectVel(dt) {
        // Semi-Lagrangian advection of velocity field
        this.newU.set(this.u);
        this.newV.set(this.v);
        const n = this.numY;
        const h = this.h;

        for (let i = 1; i < this.numX; i++) {
            for (let j = 1; j < this.numY; j++) {
                // Advect horizontal velocity
                if (this.s[i * n + j] != 0.0 && this.s[(i - 1) * n + j] != 0.0 && j < this.numY - 1) {
                    let x = i * h;
                    let y = j * h + 0.5 * h;
                    let u = this.u[i * n + j];
                    let v = (this.v[(i - 1) * n + j] + this.v[i * n + j] + this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]) * 0.25;

                    x = x - dt * u;
                    y = y - dt * v;
                    u = this.fieldCalc(x, y, U_FIELD);

                    this.newU[i * n + j] = u;
                }
                // Advect vertical velocity
                if (this.s[i * n + j] != 0.0 && this.s[i * n + j - 1] != 0.0 && i < this.numX - 1) {
                    let x = i * h + 0.5 * h;
                    let y = j * h;
                    let u = (this.u[i * n + j - 1] + this.u[i * n + j] + this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]) * 0.25;

                    let v = this.v[i * n + j];
                    x = x - dt * u;
                    y = y - dt * v;
                    v = this.fieldCalc(x, y, V_FIELD);

                    this.newV[i * n + j] = v;
                }
            }
        }
        this.u.set(this.newU);
        this.v.set(this.newV);
    }

    advectSmoke(dt) {
        // Semi-Lagrangian advection of smoke density
        this.newSMO.set(this.smo);
        const n = this.numY;
        const h = this.h;

        for (let i = 1; i < this.numX - 1; i++) {
            for (let j = 1; j < this.numY - 1; j++) {
                if (this.s[i * n + j] != 0.0) {
                    const u = (this.u[i * n + j] + this.u[(i + 1) * n + j]) * 0.5;
                    const v = (this.v[i * n + j] + this.v[i * n + j + 1]) * 0.5;
                    const x = i * h + 0.5 * h - dt * u;
                    const y = j * h + 0.5 * h - dt * v;
                    this.newSMO[i * n + j] = this.fieldCalc(x, y, S_FIELD);
                }
            }
        }
        this.smo.set(this.newSMO);
    }

    main(dt, maxIters) {
        this.p.fill(0.0); // Reset pressure
        this.solveIncomp(dt, maxIters);
        this.applyBoundaryConditions();
        this.advectVel(dt);
        this.advectSmoke(dt);
    }
}

const scene = {
    dt: 1.0 / 60.0,
    maxIters: 50,
    overRelaxation: 1.75,
    obstacleRadius: 0.1,
    obstacleX: 0.6,
    obstacleY: 0.05,
    paused: false,
    showStreamlines: false,
    showPressure: false,
    showSmoke: true,
    hiRes: false,
};

// -------------------------Canvas Creation-----------------------------
function setupScene() {
    // Initialize simulation parameters and canvas settings
    let res;
    if (scene.hiRes) {
        res = 240;
        scene.dt = 1.0 / 180.0;
    } else {
        res = 80;
    }
    const domainHeight = 1.0;
    const domainWidth = domainHeight / simHeight * simWidth;
    const h = domainHeight / res;
    const numX = Math.floor(domainWidth / h);
    const numY = Math.floor(domainHeight / h);
    const f = scene.fluid = new Fluid(numX, numY, h);
    const n = f.numY;

    // Set up inflow winds
    scene.windVel = document.getElementById('windValue').textContent;
    for (let i = 0; i < f.numX; i++) {
        for (let j = 0; j < f.numY; j++) {
            let s = 1.0;
            if (i == 0 || j == 0 || j == f.numY - 1)
                s = 0.0;
            f.s[i * n + j] = s;
            if (i == 1) {
                f.u[i * n + j] = scene.windVel;
            }
        }
    }

    // Set up smoke strip
    const smokeWidth = 0.18 * f.numY;
    const smokeOffset = -0.35;
    const minJ = Math.floor((0.5 + smokeOffset) * f.numY - 0.5 * smokeWidth);
    const maxJ = Math.floor((0.5 + smokeOffset) * f.numY + 0.5 * smokeWidth);
    for (let j = minJ; j < maxJ; j++)
        f.smo[j] = 0.0;

    setObstacle(scene.obstacleX, scene.obstacleY, true);

    // UI sync
    document.getElementById("streamButton").checked = scene.showStreamlines;
    document.getElementById("pressureButton").checked = scene.showPressure;
    document.getElementById("smokeButton").checked = scene.showSmoke;
    document.getElementById("overrelaxButton").checked = scene.overRelaxation;
    document.getElementById("hiResButton").checked = scene.hiRes;
}

// -------------------------Color Style-----------------------------
function getSciColor(val, minVal, maxVal) {
    // Normalize field values to a range of [0,1]
    val = Math.min(Math.max(val, minVal), maxVal - Number.EPSILON);
    let d = maxVal - minVal;
    val = d == 0.0 ? 0.5 : (val - minVal) / d;

    // Color into 4 segments
    const segmentSize = 0.25;
    const segmentIndex = Math.floor(val / segmentSize);
    const segmentPosition = (val - segmentIndex * segmentSize) / segmentSize;
    let color = {
        r: 0,
        g: 0,
        b: 0
    };
    switch (segmentIndex) {
        case 0:
            color.r = 0.0;
            color.g = segmentPosition;
            color.b = 1.0;
            break;
        case 1:
            color.r = 0.0;
            color.g = 1.0;
            color.b = 1.0 - segmentPosition;
            break;
        case 2:
            color.r = segmentPosition;
            color.g = 1.0;
            color.b = 0.0;
            break;
        case 3:
            color.r = 1.0;
            color.g = 1.0 - segmentPosition;
            color.b = 0.0;
            break;
    }
    return [255 * color.r, 255 * color.g, 255 * color.b, 255]
}

// ----------------- Coordinate Conversion ---------------------------
function cX(x) {
    return x * cScale;
}

function cY(y) {
    return canvas.height - y * cScale;
}

// ----------------- Map Sims Data onto Canvas ------------------------
function map() {
    // Clear the canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#DDDDDD";
    f = scene.fluid;
    const n = f.numY;
    const h = f.h;
    const cellScale = scene.hiRes ? 3.3 : 1.1; // Adjust cell scale based on resolution

    // Calculate pressure range
    minP = f.p[0];
    maxP = f.p[0];
    for (let i = 0; i < f.numCells; i++) {
        minP = Math.min(minP, f.p[i]);
        maxP = Math.max(maxP, f.p[i]);
    }

    // Get image data for pixel manipulation
    imgdata = ctx.getImageData(0, 0, canvas.width, canvas.height)
    let color = [255, 255, 255, 255]

    // Main loop for drawing fluid cells
    for (let i = 0; i < f.numX; i++) {
        for (let j = 0; j < f.numY; j++) {
            // Determine cell color based on visualization mode
            if (scene.showPressure) {
                const p = f.p[i * n + j];
                const s = f.smo[i * n + j];
                color = getSciColor(p, minP, maxP);
                if (scene.showSmoke) {
                    // Blend pressure color with smoke
                    color[0] = Math.max(0.0, color[0] - 255 * s);
                    color[1] = Math.max(0.0, color[1] - 255 * s);
                    color[2] = Math.max(0.0, color[2] - 255 * s);
                }
            } else if (scene.showSmoke) {
                // Show smoke density
                const s = f.smo[i * n + j];
                color[0] = color[1] = color[2] = 255 * s;
            } else if (f.s[i * n + j] == 0.0) {
                color[0] = color[1] = color[2] = 0; // Color for solid cells
            }

            // Calculate pixel position and size for the cell
            const x = Math.floor(cX(i * h));
            const y = Math.floor(cY((j + 1) * h));
            const cx = Math.floor(cScale * cellScale * h) + 1;
            const cy = Math.floor(cScale * cellScale * h) + 1;

            // Draw the cell pixel by pixel
            for (let yi = y; yi < y + cy; yi++) {
                let pixel = 4 * (yi * canvas.width + x)
                for (let xi = 0; xi < cx; xi++) {
                    imgdata.data[pixel++] = color[0];
                    imgdata.data[pixel++] = color[1];
                    imgdata.data[pixel++] = color[2];
                    imgdata.data[pixel++] = 255; // Alpha Value
                }
            }
        }
    }
    ctx.putImageData(imgdata, 0, 0); // Apply the pixel data to the canvas

    // Draw streamlines if enabled
    if (scene.showStreamlines) {
        const numSegs = 15;
        const lenSegs = 0.0075;
        ctx.strokeStyle = "#000000";
        const streamlineStep = scene.hiRes ? 15 : 5; // Adjust streamline density based on resolution

        // Loop through grid to draw streamlines
        for (let i = 1; i < f.numX - 1; i += streamlineStep) {
            for (let j = 1; j < f.numY - 1; j += streamlineStep) {
                let x = (i + 0.5) * f.h;
                let y = (j + 0.5) * f.h;

                ctx.beginPath();
                ctx.moveTo(cX(x), cY(y));

                // Draw individual streamline
                for (let n = 0; n < numSegs; n++) {
                    let u = f.fieldCalc(x, y, U_FIELD);
                    let v = f.fieldCalc(x, y, V_FIELD);
                    l = Math.sqrt(u * u + v * v);
                    x += u * lenSegs;
                    y += v * lenSegs;
                    if (x > f.numX * f.h)
                        break;

                    ctx.lineTo(cX(x), cY(y));
                }
                ctx.stroke();
            }
        }
    }

    // Display pressure range if pressure visualization is enabled
    if (scene.showPressure) {
        const s = "Pressure: " + minP.toFixed(1) + " - " + maxP.toFixed(1) + " Pa";
        ctx.fillStyle = scene.showSmoke ? "#FFFFFF" : "#000000";
        ctx.font = "16px Times New Roman";
        ctx.fillText(s, 10, 35);
    }
}

function drawObstacle() {
    let bottomCorrection = scene.hiRes ? 1 : 0.2;
    r = scene.obstacleRadius + f.h;

    ctx.fillStyle = scene.showPressure ? "#000000" : "#BBBBBB";
    ctx.beginPath();
    ctx.arc(cX(scene.obstacleX), cY(scene.obstacleY), cScale * r, Math.PI, 2.0 * Math.PI);
    ctx.closePath();
    ctx.fill();

    ctx.lineWidth = 4.0;
    ctx.strokeStyle = "#000000";
    ctx.beginPath();
    ctx.arc(cX(scene.obstacleX), cY(scene.obstacleY), cScale * r, Math.PI, 2.0 * Math.PI);
    ctx.closePath();
    ctx.stroke();

    ctx.lineWidth = 8.0;
    ctx.beginPath();
    ctx.moveTo(cX(scene.obstacleX - r), cY(scene.obstacleY + bottomCorrection * f.h));
    ctx.lineTo(cX(scene.obstacleX + r), cY(scene.obstacleY + bottomCorrection * f.h));
    ctx.stroke();
    ctx.lineWidth = 1.0;
}

function setObstacle(x, y, reset) {
    let vx = 0.0;
    let vy = 0.0;

    if (!reset) {
        vx = (x - scene.obstacleX) / scene.dt;
        vy = (y - scene.obstacleY) / scene.dt;
    }
    scene.obstacleX = x;
    scene.obstacleY = y;
    const r = scene.obstacleRadius;
    const f = scene.fluid;
    const n = f.numY;

    for (let i = 1; i < f.numX - 2; i++) {
        for (let j = 1; j < f.numY - 2; j++) {
            f.s[i * n + j] = 1.0;
            dx = (i + 0.5) * f.h - x;
            dy = (j + 0.5) * f.h - y;

            // Half circle, downwards facing
            if (dx * dx + dy * dy < r * r && dy > 0) {
                f.s[i * n + j] = 0.0;
                f.smo[i * n + j] = 1.0;
                f.u[i * n + j] = vx;
                f.u[(i + 1) * n + j] = vx;
                f.v[i * n + j] = vy;
                f.v[i * n + j + 1] = vy;
            }
        }
    }
}

// Interactions ------------------------------------------
let mouseDown = false;

function startDrag(x, y) {
    let bounds = canvas.getBoundingClientRect();
    let mx = x - bounds.left - canvas.clientLeft;
    let my = y - bounds.top - canvas.clientTop;
    mouseDown = true;
    x = mx / cScale;
    y = (canvas.height - my) / cScale;
    setObstacle(x, y, true);
}

function drag(x, y) {
    if (mouseDown) {
        let bounds = canvas.getBoundingClientRect();
        let mx = x - bounds.left - canvas.clientLeft;
        let my = y - bounds.top - canvas.clientTop;
        x = mx / cScale;
        y = (canvas.height - my) / cScale;
        setObstacle(x, y, false);
    }
}

function endDrag() {
    mouseDown = false;
}

canvas.addEventListener('mousedown', event => {
    startDrag(event.x, event.y);
});

canvas.addEventListener('mouseup', event => {
    endDrag();
});

canvas.addEventListener('mousemove', event => {
    drag(event.x, event.y);
});

canvas.addEventListener('touchstart', event => {
    startDrag(event.touches[0].clientX, event.touches[0].clientY)
});

canvas.addEventListener('touchend', event => {
    endDrag()
});

canvas.addEventListener('touchmove', event => {
    event.preventDefault();
    event.stopImmediatePropagation();
    drag(event.touches[0].clientX, event.touches[0].clientY)
}, {
    passive: false
});

function togglePause() {
    let button = document.getElementById('pauseButton');
    if (scene.paused)
        button.innerHTML = "Pause";
    else
        button.innerHTML = "Resume";
    scene.paused = !scene.paused;
}

function stepMotion() {
    scene.paused = false;
    simulate();
    scene.paused = true;
}

//--------------Main Simulation Loop--------------------------------
function simulate() {
    if (!scene.paused)
        scene.fluid.main(scene.dt, scene.maxIters)
}

function update() {
    simulate();
    map();
    drawObstacle();
    requestAnimationFrame(update);
}

function init() {
    setupScene();
    try {
        // Radius slider setup
        scene.obstacleRadius = parseFloat(radiusSlider.value);
        radiusValue.textContent = scene.obstacleRadius;

        radiusSlider.addEventListener('input', function() {
            scene.obstacleRadius = parseFloat(this.value);
            radiusValue.textContent = this.value;
            setObstacle(scene.obstacleX, scene.obstacleY, true);
        });

        // Wind slider setup
        windSlider.addEventListener('input', function() {
            scene.windVel = parseFloat(windSlider.value);
            windValue.textContent = scene.windVel;

            // Update the fluid simulation with new wind velocity
            for (let j = 0; j < scene.fluid.numY; j++) {
                scene.fluid.u[scene.fluid.numY + j] = scene.windVel;
            }
        });

        // HiRes checkbox setup
        let hiResButton = document.getElementById('hiResButton');
        hiResButton.checked = scene.hiRes;

        hiResButton.addEventListener('change', function() {
            scene.hiRes = this.checked;
            setupScene();
        });

        update(); // Initial update to ensure everything is drawn

    } catch (error) {
        console.error("Error setting up UI controls:", error);
    }
}
//------------- Action on Page Loading -----------
init();
