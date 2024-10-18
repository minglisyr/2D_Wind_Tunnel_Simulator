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
    {x: -0.10435750764790764, y: -0.0007323943661971832},
    {x: -0.10100655021834062, y: 0.0030140845070422535},
    {x: -0.08148999454191033, y: 0.010807453416149069},
    {x: -0.06767411428571429, y: 0.014600938967136149},
    {x: -0.0529101076426265, y: 0.015948356807511736},
    {x: -0.038234063260340634, y: 0.016600938967136148},
    {x: -0.028085424659864995, y: 0.020133802816901408},
    {x: -0.008312103174603174, y: 0.030234741784037555},
    {x: 0.0021424659864996545, y: 0.033276056338028175},
    {x: 0.016737229437229437, y: 0.034324882629107984},
    {x: 0.034945424659864995, y: 0.032535211267605634},
    {x: 0.04892465986499655, y: 0.028356807511737087},
    {x: 0.06280465986499654, y: 0.023610328638497652},
    {x: 0.08075210317460317, y: 0.020824414715719064},
    {x: 0.0922941076426265, y: 0.02156807511737089},
    {x: 0.08847229437229438, y: 0.018432863849765257},
    {x: 0.08924659864996548, y: 0.014075117370892018},
    {x: 0.09125502183406113, y: 0.007063849765258216},
    {x: 0.09147094736842105, y: -0.00029673590504451037},
    {x: 0.09499454191033139, y: -0.002996713615023474},
    {x: 0.09103722943722943, y: -0.012319248826291079},
    {x: 0.08655076479076479, y: -0.014887323943661972},
    {x: 0.0719941076426265, y: -0.016758217270194986},
    {x: 0.06341210317460317, y: -0.023334272300469483},
    {x: 0.05784465986499654, y: -0.024820610328638497},
    {x: 0.04979076479076479, y: -0.023899530516431924},
    {x: 0.041469783549783546, y: -0.019504694835680752},
    {x: 0.022963376623376622, y: -0.019375586854460093},
    {x: -0.05421502183406113, y: -0.020333802816901408},
    {x: -0.05948999454191033, y: -0.023816901408450704},
    {x: -0.06675756613756614, y: -0.025123004694835682},
    {x: -0.07351661375661375, y: -0.024123004694835682},
    {x: -0.08000655021834062, y: -0.020244131455399062},
    {x: -0.10178085424659866, y: -0.019288732394366196},
    {x: -0.1063941076426265, y: -0.017157746478873238},
    {x: -0.10378085424659866, y: -0.012708920187793428},
    {x: -0.1049551583011583, y: -0.003258215962441315}
];

const carShape_F150 = [
    {x: -0.09607843137254901, y: -0.013693438320209974},
    {x: -0.09637254901960785, y: -0.011229508196721312},
    {x: -0.09598039215686274, y: -0.008767123287671232},
    {x: -0.09499999999999999, y: -0.007067138096411856},
    {x: -0.09333333333333334, y: -0.005399543851732738},
    {x: -0.09303921568627451, y: -0.004122950819672131},
    {x: -0.0938235294117647, y: 0.0004908600228832952},
    {x: -0.09372549019607843, y: 0.0030431266846361185},
    {x: -0.09303921568627451, y: 0.005202733485193622},
    {x: -0.09176470588235295, y: 0.007067138096411856},
    {x: -0.09, y: 0.008343729128735634},
    {x: -0.08686274509803922, y: 0.009424083769633508},
    {x: -0.08215686274509804, y: 0.010406695226438188},
    {x: -0.07627450980392156, y: 0.011091958544546726},
    {x: -0.06264705882352941, y: 0.012172313185445801},
    {x: -0.05715686274509804, y: 0.013154924642249282},
    {x: -0.05284313725490196, y: 0.014627239143729128},
    {x: -0.048921568627450984, y: 0.016591186012515644},
    {x: -0.036078431372549024, y: 0.024835130970724192},
    {x: -0.030098039215686274, y: 0.027997084548104956},
    {x: -0.025, y: 0.029962486598442096},
    {x: -0.021078431372549018, y: 0.030945098055246776},
    {x: -0.016274509803921568, y: 0.03143581426611808},
    {x: 0.0025490196078431372, y: 0.03192653047698939},
    {x: 0.017450980392156862, y: 0.03231870206898939},
    {x: 0.024803921568627452, y: 0.03192653047698939},
    {x: 0.028725490196078432, y: 0.031239811777751556},
    {x: 0.031274509803921568, y: 0.030357744297384184},
    {x: 0.032745098039215686, y: 0.029178846239713724},
    {x: 0.03323529411764706, y: 0.027997084548104956},
    {x: 0.033725490196078434, y: 0.017277904711753478},
    {x: 0.034313725490196076, y: 0.013940083313905872},
    {x: 0.035098039215686274, y: 0.012663492281582094},
    {x: 0.036470588235294116, y: 0.0120750668692038},
    {x: 0.040392156862745096, y: 0.011977978087874508},
    {x: 0.053431372549019608, y: 0.012762036444040678},
    {x: 0.07264705882352942, y: 0.01305566541736997},
    {x: 0.09009803921568627, y: 0.013745295333937096},
    {x: 0.09441176470588236, y: 0.013451666360607804},
    {x: 0.09637254901960785, y: 0.012861494566499204},
    {x: 0.09745098039215687, y: 0.011977978087874508},
    {x: 0.09803921568627451, y: 0.01050420620506537},
    {x: 0.09794117647058824, y: 0.006870960533755272},
    {x: 0.09686274509803922, y: -0.00373029358600583},
    {x: 0.09696078431372549, y: -0.007067138096411856},
    {x: 0.09745098039215687, y: -0.008638579125012648},
    {x: 0.09813725490196078, y: -0.009030750717012648},
    {x: 0.09931372549019608, y: -0.009226928279669232},
    {x: 0.09980392156862746, y: -0.010307282920567108},
    {x: 0.09990196078431372, y: -0.012861494566499204},
    {x: 0.09941176470588236, y: -0.014530544400507616},
    {x: 0.09862745098039215, y: -0.015513155857312296},
    {x: 0.09715686274509804, y: -0.016415612337679668},
    {x: 0.09431372549019608, y: -0.017005784131788268},
    {x: 0.08107843137254901, y: -0.017988395588592948},
    {x: 0.07499999999999999, y: -0.017790760644607072},
    {x: 0.07343137254901961, y: -0.018282932236807072},
    {x: 0.07254901960784313, y: -0.019265543693611752},
    {x: 0.07176470588235294, y: -0.021917391756887232},
    {x: 0.07098039215686274, y: -0.024851517125542712},
    {x: 0.06970588235294117, y: -0.026913640413223236},
    {x: 0.06794117647058824, y: -0.028468318496030376},
    {x: 0.06549019607843137, y: -0.030038630997508924},
    {x: 0.06284313725490196, y: -0.031021242454313604},
    {x: 0.06, y: -0.031611414248422204},
    {x: 0.05715686274509804, y: -0.031611414248422204},
    {x: 0.05441176470588236, y: -0.030856329967855076},
    {x: 0.05186274509803922, y: -0.029579738935531298},
    {x: 0.048823529411764706, y: -0.027109817975879822},
    {x: 0.044313725490196076, y: -0.02297393698418197},
    {x: 0.042254901960784314, y: -0.02179503892651151},
    {x: 0.039705882352941174, y: -0.02120486713240291},
    {x: 0.03470588235294118, y: -0.02071269554020291},
    {x: 0.015490196078431373, y: -0.02120486713240291},
    {x: -0.006764705882352941, y: -0.021107778351073618},
    {x: -0.034705882352941176, y: -0.02179503892651151},
    {x: -0.048431372549019606, y: -0.021598861363854926},
    {x: -0.05109803921568628, y: -0.022091032956054926},
    {x: -0.05235294117647058, y: -0.022778293531492818},
    {x: -0.0538235294117647, y: -0.024638333236054926},
    {x: -0.05578431372549018, y: -0.027011273813221238},
    {x: -0.05862745098039215, y: -0.029384214390387548},
    {x: -0.06156862745098039, y: -0.031052875386524368},
    {x: -0.06382352941176471, y: -0.031611414248422204},
    {x: -0.06764705882352942, y: -0.031708503029751496},
    {x: -0.07098039215686274, y: -0.031316331437751496},
    {x: -0.07313725490196079, y: -0.030431359577797508},
    {x: -0.07549019607843137, y: -0.028762309743789096},
    {x: -0.07745098039215687, y: -0.026617646839893946},
    {x: -0.07901960784313726, y: -0.024344704262725634},
    {x: -0.08019607843137254, y: -0.022091032956054926},
    {x: -0.08107843137254902, y: -0.021715950146725634},
    {x: -0.08305882352941176, y: -0.02179503892651151},
    {x: -0.08647058823529411, y: -0.022581755988407930},
    {x: -0.08715686274509804, y: -0.022286671633749346},
    {x: -0.08764705882352941, y: -0.019930819839252328},
    {x: -0.08852941176470589, y: -0.019537190865923036},
    {x: -0.09147058823529411, y: -0.019438646703264452},
    {x: -0.09284313725490196, y: -0.018848474909155852},
    {x: -0.09421568627450981, y: -0.017573341258161366},
    {x: -0.09549019607843138, y: -0.015513155857312296}
];

const carShape_Minivan = [
    {x: -0.09995824634655532, y: 0.0006830699881461179},
    {x: -0.09912316715542522, y: 0.002951235596690735},
    {x: -0.09608351328307094, y: 0.00704523072683938},
    {x: -0.09512733082706766, y: 0.009473396335383996},
    {x: -0.09369937202288548, y: 0.01430772755247322},
    {x: -0.09289540284360189, y: 0.01572579316102262},
    {x: -0.09185032958801498, y: 0.01664983596529732},
    {x: -0.08925678496868476, y: 0.018190934171840324},
    {x: -0.08541962316715542, y: 0.02127313058493633},
    {x: -0.08285718854431599, y: 0.022537180791771335},
    {x: -0.07834655532359082, y: 0.02405524640032073},
    {x: -0.07007929515418502, y: 0.02609135481314483},
    {x: -0.06574634655532359, y: 0.02775744302597893},
    {x: -0.06162733082706766, y: 0.03012557404307773},
    {x: -0.05616715542521994, y: 0.033952706859365516},
    {x: -0.048477386934673365, y: 0.03991592528620811},
    {x: -0.04364357804719848, y: 0.04315110109502941},
    {x: -0.038658801498127344, y: 0.04624329750812542},
    {x: -0.033433832335329344, y: 0.04879338572094952},
    {x: -0.028125, y: 0.05082949413377362},
    {x: -0.022666144200626958, y: 0.05219557234660771},
    {x: -0.015274059561128526, y: 0.05340560515088241},
    {x: -0.006027855228281609, y: 0.05424767075943181},
    {x: 0.014865418502202645, y: 0.05508973636798121},
    {x: 0.031897676579925652, y: 0.05542376397225591},
    {x: 0.05166457023060797, y: 0.05491772116798121},
    {x: 0.07434048257372654, y: 0.0535516429551471},
    {x: 0.07882001879699248, y: 0.05304560015087241},
    {x: 0.07898127340823971, y: 0.05271157254659771},
    {x: 0.07737577639751553, y: 0.04966739873777641},
    {x: 0.07721452178626831, y: 0.04813933312922701},
    {x: 0.07761828120301129, y: 0.04644523011211821},
    {x: 0.07857446365901457, y: 0.04458314189928411},
    {x: 0.08351328307094673, y: 0.03991592528620811},
    {x: 0.08615794200626959, y: 0.03598879247091942},
    {x: 0.08832957679738562, y: 0.03206165965563073},
    {x: 0.08944701003135857, y: 0.02876850644108413},
    {x: 0.09017193370165746, y: 0.02547535322653753},
    {x: 0.09033318831290468, y: 0.02272526501371343},
    {x: 0.08952921913362109, y: 0.01407171755247322},
    {x: 0.08969047374486831, y: 0.01094551893937722},
    {x: 0.09017193370165746, y: 0.00894142532082782},
    {x: 0.09177743071238163, y: 0.00536626990126442},
    {x: 0.09217993994200626, y: 0.003285263201005334},
    {x: 0.09209993994200626, y: 1.4997179291917992e-5},
    {x: 0.09105486668641935, y: -0.005947255444265283},
    {x: 0.09009868423041607, y: -0.00940837626735708},
    {x: 0.08914250177441279, y: -0.011104457685926478},
    {x: 0.08786868423041607, y: -0.011939547093495277},
    {x: 0.08627444744200626, y: -0.012357591501064876},
    {x: 0.08158445072816715, y: -0.012775635908634475},
    {x: 0.07714326205450733, y: -0.014028769919487673},
    {x: 0.07361204470852018, y: -0.014863859327056472},
    {x: 0.0706308288288288, y: -0.014946864930341272},
    {x: 0.06776961295913682, y: -0.014863859327056472},
    {x: 0.06672453970354991, y: -0.01528190373462607},
    {x: 0.06600087625556738, y: -0.01611799314219487},
    {x: 0.06487218279928078, y: -0.021245156362466871},
    {x: 0.06391600034327751, y: -0.02299829037333567},
    {x: 0.06239920287259951, y: -0.02491835158841927},
    {x: 0.06004886854431599, y: -0.027271485599217672},
    {x: 0.05745532392498576, y: -0.029040647205730672},
    {x: 0.05474047884569169, y: -0.02980767960798387},
    {x: 0.05190563376639763, y: -0.02980767960798387},
    {x: 0.04927223936937577, y: -0.029040647205730672},
    {x: 0.045741022023387616, y: -0.027537531103402872},
    {x: 0.04244980515103333, y: -0.025453314392764072},
    {x: 0.04067175268035533, y: -0.023535997584525472},
    {x: 0.039635679424768426, y: -0.021701791366102673},
    {x: 0.03843829596848183, y: -0.018616523944394875},
    {x: 0.037548613512478554, y: -0.017781434536826076},
    {x: 0.03634997983387832, y: -0.017446406932511477},
    {x: 0.029666144200626958, y: -0.017363401329226677},
    {x: 0.006833832335329341, y: -0.017781434536826076},
    {x: -0.043402178626831735, y: -0.017446406932511477},
    {x: -0.05359991715542522, y: -0.017864451340081076},
    {x: -0.05497495072816715, y: -0.018199478944395676},
    {x: -0.05561991715542522, y: -0.018784517748680475},
    {x: -0.05626488358268329, y: -0.020827111954897272},
    {x: -0.05795982035336129, y: -0.02324930157990527},
    {x: -0.060697165432655355, y: -0.02491835158841927},
    {x: -0.06303624953862523, y: -0.027271485599217672},
    {x: -0.06559867438026831, y: -0.029040647205730672},
    {x: -0.06831351945956237, y: -0.02980767960798387},
    {x: -0.07114836453885644, y: -0.02997369081455347},
    {x: -0.07398320961815051, y: -0.02964167400469917},
    {x: -0.07637465371681039, y: -0.028791635999045872},
    {x: -0.07813144596517475, y: -0.027537531103402872},
    {x: -0.07956066499166957, y: -0.02580233680604847},
    {x: -0.08142746724003393, y: -0.021245156362466871},
    {x: -0.08246128027330719, y: -0.019574674737826873},
    {x: -0.08341746272931047, y: -0.018950529558535075},
    {x: -0.08461484640791069, y: -0.018784517748680475},
    {x: -0.08817741875390154, y: -0.018867523351965275},
    {x: -0.09265794200626959, y: -0.018199478944395676},
    {x: -0.09580885610093752, y: -0.017781434536826076},
    {x: -0.09701749999999999, y: -0.016946345129257277},
    {x: -0.09789498391596711, y: -0.015778965537880269},
    {x: -0.09877246783193423, y: -0.013861758109633073},
    {x: -0.09932869448150594, y: -0.011522502685925678},
    {x: -0.10029613715982286, y: -0.002444932426842685},
    {x: -0.10013488254857564, y: -0.00044304722827768}
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
