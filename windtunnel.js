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
    overRelaxation: 1.95,
    obstacleRadius: 0.1,
    obstacleX: 0.6,
    obstacleY: 0.1,
    paused: false,
    showStreamlines: false,
    showPressure: false,
    showSmoke: true,
    hiRes: false,
};

// Ferrari carshape
const carShape_Ferrari = [
    {x: -0.5217875382395382, y: -0.003661971830985916},
    {x: -0.5050327510917031, y: 0.015070422535211268},
    {x: -0.40744997270955165, y: 0.054037267080745345},
    {x: -0.33837057142857145, y: 0.07300469483568075},
    {x: -0.2645505382131325, y: 0.07974178403755868},
    {x: -0.19117031630170316, y: 0.08300469483568075},
    {x: -0.14042712329932498, y: 0.10066901408450704},
    {x: -0.04156051587301587, y: 0.15117370892018778},
    {x: 0.010712329932498273, y: 0.16638028169014087},
    {x: 0.08368614718614719, y: 0.17162441314553992},
    {x: 0.17472712329932498, y: 0.16267605633802816},
    {x: 0.24462329932498274, y: 0.14178403755868544},
    {x: 0.3140232993249827, y: 0.11805164319248826},
    {x: 0.40376051587301585, y: 0.10412207357859532},
    {x: 0.4614705382131325, y: 0.10784037558685446},
    {x: 0.4423614718614719, y: 0.09216431924882629},
    {x: 0.4462329932498274, y: 0.07037558685446009},
    {x: 0.45627510917030566, y: 0.03531924882629108},
    {x: 0.45735473684210525, y: -0.0014836795252225519},
    {x: 0.47497270955165693, y: -0.014983568075117372},
    {x: 0.45518614718614717, y: -0.061596244131455396},
    {x: 0.43275382395382395, y: -0.07443661971830986},
    {x: 0.3599705382131325, y: -0.08379108635097493},
    {x: 0.31706051587301586, y: -0.11667136150234742},
    {x: 0.2892232993249827, y: -0.12410305164319249},
    {x: 0.24895382395382395, y: -0.11949765258215962},
    {x: 0.20734891774891774, y: -0.09752347417840376},
    {x: 0.11481688311688311, y: -0.09687793427230047},
    {x: -0.27107510917030566, y: -0.10166901408450704},
    {x: -0.29744997270955164, y: -0.11908450704225352},
    {x: -0.3337878306878307, y: -0.12561502347417841},
    {x: -0.36758306878306877, y: -0.12061502347417841},
    {x: -0.4000327510917031, y: -0.10122065727699531},
    {x: -0.5089042712329933, y: -0.09644366197183099},
    {x: -0.5319705382131325, y: -0.08578873239436619},
    {x: -0.5189042712329933, y: -0.06354460093896714},
    {x: -0.5247757915057915, y: -0.016291079812206575}
];

// Ford F150 carshape
const carShape_F150 = [
    {x: -0.4804901960784314, y: -0.21527777777777777},
    {x: -0.4735294117647059, y: -0.10524691358024691},
    {x: -0.4661764705882353, y: -0.04953703703703704},
    {x: -0.4642156862745098, y: 0.08825617283950617},
    {x: -0.4392156862745098, y: 0.14399691358024691},
    {x: -0.2686274509803922, y: 0.22453703703703703},
    {x: -0.12107843137254903, y: 0.47685185185185186},
    {x: 0.11813725490196078, y: 0.5046296296296297},
    {x: 0.16421568627450982, y: 0.45679012345679015},
    {x: 0.17647058823529413, y: 0.19675925925925927},
    {x: 0.48186274509803924, y: 0.20293209876543212},
    {x: 0.49019607843137253, y: 0.13475308641975309},
    {x: 0.4892156862745098, y: -0.14089506172839506},
    {x: 0.4995098039215686, y: -0.17808641975308643},
    {x: 0.48774509803921567, y: -0.25555555555555554},
    {x: 0.3637254901960784, y: -0.2987654320987654},
    {x: 0.34313725490196076, y: -0.44},
    {x: 0.3, y: -0.4987654320987654},
    {x: 0.2554901960784314, y: -0.45679012345679015},
    {x: 0.20735294117647058, y: -0.34074074074074073},
    {x: -0.25980392156862747, y: -0.35462962962962963},
    {x: -0.31029411764705884, y: -0.49259259259259256},
    {x: -0.36568627450980393, y: -0.48024691358024694},
    {x: -0.4102941176470588, y: -0.34074074074074073},
    {x: -0.43676470588235294, y: -0.34537037037037035},
    {x: -0.4460784313725491, y: -0.30662037037037035},
    {x: -0.47254901960784315, y: -0.2709876543209877}
];

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
const carShape = carShape_F150;

// Function to draw the car
function drawObstacle(centerX, centerY, scaleRadius) {
    ctx.lineWidth = 4.0;
    ctx.beginPath();
    ctx.moveTo(
        cX(centerX + carShape[0].x * scaleRadius),
        cY(centerY + carShape[0].y * scaleRadius)
    );

    for (let i = 1; i < carShape.length; i++) {
        ctx.lineTo(
            cX(centerX + carShape[i].x * scaleRadius),
            cY(centerY + carShape[i].y * scaleRadius)
        );
    }

    ctx.closePath();
    ctx.strokeStyle = '#000000';
    ctx.fillStyle = scene.showPressure ? "#000000" : "#FF0000";
    ctx.fill();
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
    const f = scene.fluid;
    const n = f.numY;

    // Calculate bounding box of the car shape
    let minX = Infinity,
        maxX = -Infinity,
        minY = Infinity,
        maxY = -Infinity;
    for (let point of carShape) {
        minX = Math.min(minX, point.x);
        maxX = Math.max(maxX, point.x);
        minY = Math.min(minY, point.y);
        maxY = Math.max(maxY, point.y);
    }

    // Scale factor (adjust this to change the size of the car)
    const scale = scene.obstacleRadius * 2 / Math.max(maxX - minX, maxY - minY);

    for (let i = 1; i < f.numX - 2; i++) {
        for (let j = 1; j < f.numY - 2; j++) {
            f.s[i * n + j] = 1.0;
            let dx = (i + 0.5) * f.h - x;
            let dy = (j + 0.5) * f.h - y;

            // Transform to car's local coordinates
            let localX = dx / scale;
            let localY = dy / scale;

            if (isPointInCarShape(localX, localY)) {
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

function isPointInCarShape(x, y) {
    let inside = false;
    for (let i = 0, j = carShape.length - 1; i < carShape.length; j = i++) {
        let xi = carShape[i].x,
            yi = carShape[i].y;
        let xj = carShape[j].x,
            yj = carShape[j].y;

        let intersect = ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    return inside;
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
    drawObstacle(scene.obstacleX, scene.obstacleY, scene.obstacleRadius * 2.1);
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
