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

// carShapes Library
const carShape_Ferrari = [
    {x: -0.104358, y: -0.000732},
    {x: -0.101007, y: 0.003014},
    {x: -0.081490, y: 0.010807},
    {x: -0.067674, y: 0.014601},
    {x: -0.052910, y: 0.015948},
    {x: -0.038234, y: 0.016601},
    {x: -0.028085, y: 0.020134},
    {x: -0.008312, y: 0.030235},
    {x: 0.002142, y: 0.033276},
    {x: 0.016737, y: 0.034325},
    {x: 0.034945, y: 0.032535},
    {x: 0.048925, y: 0.028357},
    {x: 0.062805, y: 0.023610},
    {x: 0.080752, y: 0.020824},
    {x: 0.092294, y: 0.021568},
    {x: 0.088472, y: 0.018433},
    {x: 0.089247, y: 0.014075},
    {x: 0.091255, y: 0.007064},
    {x: 0.091471, y: -0.000297},
    {x: 0.094995, y: -0.002997},
    {x: 0.091037, y: -0.012319},
    {x: 0.086551, y: -0.014887},
    {x: 0.071994, y: -0.016758},
    {x: 0.063412, y: -0.023334},
    {x: 0.057845, y: -0.024821},
    {x: 0.049791, y: -0.023900},
    {x: 0.041470, y: -0.019505},
    {x: 0.022963, y: -0.019376},
    {x: -0.054215, y: -0.020334},
    {x: -0.059490, y: -0.023817},
    {x: -0.066758, y: -0.025123},
    {x: -0.073517, y: -0.024123},
    {x: -0.080007, y: -0.020244},
    {x: -0.101781, y: -0.019289},
    {x: -0.106394, y: -0.017158},
    {x: -0.103781, y: -0.012709},
    {x: -0.104955, y: -0.003258}
];

const carShape_F150 = [
    {x: -0.124983, y: -0.017813},
    {x: -0.125367, y: -0.014607},
    {x: -0.124867, y: -0.011407},
    {x: -0.123583, y: -0.009193},
    {x: -0.121417, y: -0.007023},
    {x: -0.121033, y: -0.005363},
    {x: -0.122050, y: 0.000638},
    {x: -0.121917, y: 0.003958},
    {x: -0.121033, y: 0.006768},
    {x: -0.119367, y: 0.009193},
    {x: -0.117083, y: 0.010855},
    {x: -0.113000, y: 0.012258},
    {x: -0.106867, y: 0.013537},
    {x: -0.099217, y: 0.014428},
    {x: -0.081483, y: 0.015833},
    {x: -0.074350, y: 0.017112},
    {x: -0.068733, y: 0.019028},
    {x: -0.063633, y: 0.021583},
    {x: -0.046933, y: 0.032305},
    {x: -0.039150, y: 0.036418},
    {x: -0.032525, y: 0.038975},
    {x: -0.027417, y: 0.040253},
    {x: -0.021167, y: 0.040892},
    {x: 0.003316, y: 0.041532},
    {x: 0.022698, y: 0.042043},
    {x: 0.032263, y: 0.041532},
    {x: 0.037367, y: 0.040637},
    {x: 0.040683, y: 0.039485},
    {x: 0.042592, y: 0.037955},
    {x: 0.043228, y: 0.036418},
    {x: 0.043867, y: 0.022478},
    {x: 0.044633, y: 0.018135},
    {x: 0.045650, y: 0.016477},
    {x: 0.047433, y: 0.015708},
    {x: 0.052537, y: 0.015582},
    {x: 0.069508, y: 0.016603},
    {x: 0.094492, y: 0.016985},
    {x: 0.117200, y: 0.017875},
    {x: 0.122808, y: 0.017493},
    {x: 0.125367, y: 0.016732},
    {x: 0.126758, y: 0.015582},
    {x: 0.127517, y: 0.013663},
    {x: 0.127400, y: 0.008935},
    {x: 0.125992, y: -0.004853},
    {x: 0.126117, y: -0.009193},
    {x: 0.126758, y: -0.011235},
    {x: 0.127658, y: -0.011745},
    {x: 0.129183, y: -0.012000},
    {x: 0.129833, y: -0.013407},
    {x: 0.129958, y: -0.016732},
    {x: 0.129317, y: -0.018898},
    {x: 0.128292, y: -0.020175},
    {x: 0.126383, y: -0.021353},
    {x: 0.122683, y: -0.022122},
    {x: 0.105458, y: -0.023398},
    {x: 0.097558, y: -0.023143},
    {x: 0.095517, y: -0.023782},
    {x: 0.094367, y: -0.025058},
    {x: 0.093350, y: -0.028508},
    {x: 0.092333, y: -0.032342},
    {x: 0.090667, y: -0.035000},
    {x: 0.088383, y: -0.037038},
    {x: 0.085183, y: -0.039077},
    {x: 0.081750, y: -0.040353},
    {x: 0.078075, y: -0.041118},
    {x: 0.074350, y: -0.041118},
    {x: 0.070775, y: -0.040125},
    {x: 0.067467, y: -0.038465},
    {x: 0.063508, y: -0.035255},
    {x: 0.057642, y: -0.029885},
    {x: 0.054967, y: -0.028353},
    {x: 0.051658, y: -0.027587},
    {x: 0.045142, y: -0.026945},
    {x: 0.020150, y: -0.027587},
    {x: -0.008798, y: -0.027458},
    {x: -0.045142, y: -0.028353},
    {x: -0.062992, y: -0.028098},
    {x: -0.066467, y: -0.028735},
    {x: -0.068100, y: -0.029628},
    {x: -0.070017, y: -0.032048},
    {x: -0.072567, y: -0.035128},
    {x: -0.076267, y: -0.038212},
    {x: -0.080083, y: -0.040382},
    {x: -0.083025, y: -0.041118},
    {x: -0.087992, y: -0.041245},
    {x: -0.092333, y: -0.040737},
    {x: -0.095133, y: -0.039587},
    {x: -0.098192, y: -0.037418},
    {x: -0.100742, y: -0.034628},
    {x: -0.102792, y: -0.031665},
    {x: -0.104317, y: -0.028735},
    {x: -0.105458, y: -0.028252},
    {x: -0.108042, y: -0.028353},
    {x: -0.112483, y: -0.029373},
    {x: -0.113367, y: -0.028992},
    {x: -0.114017, y: -0.025933},
    {x: -0.115158, y: -0.025417},
    {x: -0.118983, y: -0.025288},
    {x: -0.120767, y: -0.024523},
    {x: -0.122550, y: -0.022863},
    {x: -0.124208, y: -0.020175}
];

const carShape_Minivan = [
    {x: -0.112000, y: 0.000765},
    {x: -0.111067, y: 0.003307},
    {x: -0.107667, y: 0.007893},
    {x: -0.106600, y: 0.010093},
    {x: -0.105000, y: 0.015493},
    {x: -0.104100, y: 0.016053},
    {x: -0.102933, y: 0.016667},
    {x: -0.100000, y: 0.017053},
    {x: -0.095733, y: 0.018373},
    {x: -0.092867, y: 0.018853},
    {x: -0.087800, y: 0.019573},
    {x: -0.078533, y: 0.020000},
    {x: -0.073667, y: 0.020560},
    {x: -0.069067, y: 0.021760},
    {x: -0.062933, y: 0.023147},
    {x: -0.054333, y: 0.025573},
    {x: -0.048933, y: 0.026787},
    {x: -0.043333, y: 0.028000},
    {x: -0.037467, y: 0.028920},
    {x: -0.031533, y: 0.029720},
    {x: -0.025400, y: 0.030253},
    {x: -0.017120, y: 0.030707},
    {x: -0.006760, y: 0.031053},
    {x: 0.016653, y: 0.031333},
    {x: 0.035733, y: 0.031413},
    {x: 0.057867, y: 0.031173},
    {x: 0.083333, y: 0.030707},
    {x: 0.088333, y: 0.030573},
    {x: 0.088533, y: 0.030253},
    {x: 0.086733, y: 0.029253},
    {x: 0.086533, y: 0.028613},
    {x: 0.087000, y: 0.028000},
    {x: 0.088067, y: 0.027333},
    {x: 0.093600, y: 0.025573},
    {x: 0.096533, y: 0.024000},
    {x: 0.099000, y: 0.022427},
    {x: 0.100267, y: 0.021173},
    {x: 0.101067, y: 0.019920},
    {x: 0.101267, y: 0.018853},
    {x: 0.100333, y: 0.015493},
    {x: 0.100533, y: 0.014307},
    {x: 0.101067, y: 0.013547},
    {x: 0.102867, y: 0.012173},
    {x: 0.103333, y: 0.011373},
    {x: 0.103267, y: 0.010120},
    {x: 0.102067, y: 0.007827},
    {x: 0.101000, y: 0.006493},
    {x: 0.099933, y: 0.005840},
    {x: 0.098467, y: 0.005520},
    {x: 0.096667, y: 0.005360},
    {x: 0.091400, y: 0.005200},
    {x: 0.086467, y: 0.004720},
    {x: 0.082467, y: 0.004400},
    {x: 0.079200, y: 0.004373},
    {x: 0.075800, y: 0.004400},
    {x: 0.074733, y: 0.004240},
    {x: 0.073933, y: 0.003920},
    {x: 0.072667, y: 0.001947},
    {x: 0.071600, y: 0.001267},
    {x: 0.069933, y: 0.000533},
    {x: 0.067267, y: -0.000373},
    {x: 0.064367, y: -0.001053},
    {x: 0.061333, y: -0.001347},
    {x: 0.058167, y: -0.001347},
    {x: 0.055200, y: -0.001053},
    {x: 0.051253, y: -0.000480},
    {x: 0.047567, y: 0.000320},
    {x: 0.045567, y: 0.001067},
    {x: 0.044400, y: 0.001973},
    {x: 0.043300, y: 0.003493},
    {x: 0.042533, y: 0.004187},
    {x: 0.041467, y: 0.004507},
    {x: 0.034333, y: 0.004533},
    {x: 0.007933, y: 0.004107},
    {x: -0.050333, y: 0.004453},
    {x: -0.062200, y: 0.004080},
    {x: -0.063800, y: 0.003733},
    {x: -0.064533, y: 0.003120},
    {x: -0.065267, y: 0.000907},
    {x: -0.067267, y: -0.001627},
    {x: -0.070467, y: -0.003413},
    {x: -0.073133, y: -0.005893},
    {x: -0.076133, y: -0.007813},
    {x: -0.079267, y: -0.008627},
    {x: -0.082533, y: -0.008800},
    {x: -0.085800, y: -0.008453},
    {x: -0.088600, y: -0.007547},
    {x: -0.090600, y: -0.006213},
    {x: -0.092333, y: -0.004373},
    {x: -0.094467, y: 0.000320},
    {x: -0.095667, y: 0.002080},
    {x: -0.096733, y: 0.002747},
    {x: -0.098133, y: 0.002933},
    {x: -0.102267, y: 0.002853},
    {x: -0.107467, y: 0.003573},
    {x: -0.111133, y: 0.004027},
    {x: -0.112467, y: 0.004907},
    {x: -0.113467, y: 0.006133},
    {x: -0.114467, y: 0.008187},
    {x: -0.115133, y: 0.010613},
    {x: -0.116267, y: 0.020293},
    {x: -0.116067, y: 0.022507}
];

// setup scene
const scene = {
    dt: 1.0 / 60.0,
    maxIters: 50,
    overRelaxation: 1.95,
    carShape: carShape_Ferrari,
    obstacleRadius: 1.0,
    obstacleX: 0.6,
    obstacleY: 0.1,
    paused: false,
    showStreamlines: false,
    showPressure: false,
    showSmoke: true,
    hiRes: false,
    fluid: null
};

// switch car model
function changeCarModel(model) {
    switch (model) {
        case 'Ferrari':
            scene.carShape = carShape_Ferrari;
            break;
        case 'Ford F150':
            scene.carShape = carShape_F150;
            break;
        case 'Minivan':
            scene.carShape = carShape_Minivan;
            break;
    }
    setObstacle(scene.obstacleX, scene.obstacleY, true);
    update();
}

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

// Function to draw the car
function drawObstacle(centerX, centerY, scaleRadius) {
    ctx.lineWidth = 4.0;
    ctx.beginPath();
    ctx.moveTo(
        cX(centerX + scene.carShape[0].x * scaleRadius),
        cY(centerY + scene.carShape[0].y * scaleRadius)
    );

    for (let i = 1; i < scene.carShape.length; i++) {
        ctx.lineTo(
            cX(centerX + scene.carShape[i].x * scaleRadius),
            cY(centerY + scene.carShape[i].y * scaleRadius)
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
    for (let point of scene.carShape) {
        minX = Math.min(minX, point.x);
        maxX = Math.max(maxX, point.x);
        minY = Math.min(minY, point.y);
        maxY = Math.max(maxY, point.y);
    }

    // Scale factor (adjust this to change the size of the car)
    const scale = 0.25 * scene.obstacleRadius / Math.max(maxX - minX, maxY - minY);

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
    for (let i = 0, j = scene.carShape.length - 1; i < scene.carShape.length; j = i++) {
        let xi = scene.carShape[i].x,
            yi = scene.carShape[i].y;
        let xj = scene.carShape[j].x,
            yj = scene.carShape[j].y;

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
    drawObstacle(scene.obstacleX, scene.obstacleY, scene.obstacleRadius * 1.5);
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
document.addEventListener('DOMContentLoaded', () => {
    init();
});
