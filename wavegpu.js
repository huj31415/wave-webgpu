let adapter, device;

const ui = Object.freeze({
  panel: document.getElementById("controls"),
  collapse: document.getElementById("toggleSettings"),
  dt: document.getElementById("dt"),
  dtValue: document.getElementById("dtValue"),
  displayType: document.getElementById("visualization"),
  reInit: document.getElementById("reinit"),
  simSpeedSlider: document.getElementById("simSpeed"),
  simSpeedValue: document.getElementById("simSpeedValue"),
  imageUpload: document.getElementById("barrierUpload"),
  blackRefractSlider: document.getElementById("blackIRSlider"),
  blackRefractValue: document.getElementById("blackIRValue"),
  whiteRefractSlider: document.getElementById("whiteIRSlider"),
  whiteRefractValue: document.getElementById("whiteIRValue"),
  imageApply: document.getElementById("applyBarrierImage"),
  imageScale: document.getElementById("imageScale"),
  barrierInvert: document.getElementById("barrierInvert"),
  cClear: document.getElementById("clearBarriers"),
  boundaryAbsorb: document.getElementById("boundaryAbsorb"),
  wavelengthSlider: document.getElementById("wavelength"),
  wavelengthValue: document.getElementById("wavelengthValue"),
  amplitudeSlider: document.getElementById("amplitude"),
  amplitudeValue: document.getElementById("amplitudeValue"),
  drawnBarrierAbsorb: document.getElementById("drawnBarrierAbsorb"),
  planeWave: document.getElementById("planeWave"),
  pointSource: document.getElementById("pointSource"),
  waveToggle: document.getElementById("waveToggle"),
  blockTop: document.getElementById("blockTop"),
  blockBottom: document.getElementById("blockBottom"),
  jsTime: document.getElementById("jsTime"),
  // gpuTime: document.getElementById("gpuTime"),
  frameTime: document.getElementById("frameTime"),
  fps: document.getElementById("fps"),
  // preset settings
  preset: document.getElementById("preset"),
  presetSettings: document.getElementById("presetSettings"),
  loadPreset: document.getElementById("loadPreset"),
  offsetLeftSlider: document.getElementById("offsetLeftSlider"),
  offsetLeftValue: document.getElementById("offsetLeftValue"),
  dsSettings: document.getElementById("dsSettings"),
  dsSpacingSlider: document.getElementById("dsSpacingSlider"),
  dsSpacingValue: document.getElementById("dsSpacingValue"),
  dsWidthSlider: document.getElementById("dsWidthSlider"),
  dsWidthValue: document.getElementById("dsWidthValue"),
  dsNumSlits: document.getElementById("dsNumSlits"),
  lensSettings: document.getElementById("lensSettings"),
  lensRadiusSlider: document.getElementById("lensRadiusSlider"),
  lensRadiusValue: document.getElementById("lensRadiusValue"),
  lensARSlider: document.getElementById("lensARSlider"),
  lensARValue: document.getElementById("lensARValue"),
  lensIRSlider: document.getElementById("lensIRSlider"),
  lensIRValue: document.getElementById("lensIRValue"),
  rSettings: document.getElementById("rSettings"),
  rRadiusSlider: document.getElementById("rRadiusSlider"),
  rRadiusValue: document.getElementById("rRadiusValue"),
  rDir: document.getElementById("rDir"),
  lensOptions: document.getElementById("lensOptions"),
  lensDir: document.getElementById("lensDir"),
  lensDRSlider: document.getElementById("lensDRSlider"),
  lensDRValue: document.getElementById("lensDRValue"),
});

async function main() {

  const width = window.innerWidth;
  const height = window.innerHeight;
  const halfWidth = Math.round(width / 2);
  const halfHeight = Math.round(height / 2);

  // Simulation Parameters
  let numCells = width * height;

  let dtPerFrame = parseInt(ui.simSpeedSlider.value);
  let dt = parseFloat(ui.dt.value);//1 / Math.sqrt(2) - 1e-4; // < dx / c*sqrt2

  const initU = 0; // initial state of simulation domain
  const initC = 1; // initial value of C

  let waveOn = true;
  const decayRate = 0.1; // decay rate when wave is disabled
  let ampVal = 1;

  let exit = false;

  // WebGPU Setup
  const canvas = document.getElementById("canvas");
  canvas.width = width;
  canvas.height = height;
  if (!adapter) {
    adapter = await navigator.gpu?.requestAdapter();
    const canTimestamp = adapter.features.has('timestamp-query');
    device = await adapter?.requestDevice({
      requiredFeatures: [
        ...(canTimestamp ? ['timestamp-query'] : []),
      ],
    });
  }
  if (!device) {
    alert("Browser does not support WebGPU");
    document.body.textContent = "WebGPU is not supported in this browser.";
    return;
  }
  const context = canvas.getContext("webgpu");
  const swapChainFormat = "bgra8unorm";
  context.configure({
    device: device,
    format: swapChainFormat,
  });

  const displayType = Object.freeze({
    wave: 0,
    speed: 1,
    end: 2
  })

  // Preset setup
  const Preset = Object.freeze({
    slits: 0,
    lens: 1,
    circularReflector: 2,
    parabolicReflector: 3,
  });

  let activePreset = Preset[ui.preset.value];

  ui.offsetLeftSlider.max = width - 1;
  // ui.offsetLeftSlider.value = ui.offsetLeftValue.value = 200;
  ui.dsSpacingSlider.max = ui.dsSpacingValue.textContent = height;
  // ui.dsSpacingSlider.value = ui.dsSpacingValue.textContent = 200;
  ui.lensRadiusSlider.max = height * 2;
  // ui.lensRadiusValue.value = 300;
  ui.rRadiusSlider.max = ui.rRadiusSlider.value = ui.rRadiusValue.textContent = height;

  // Simulation Buffers
  // current, previous, and next time step
  const stateBufferSize = numCells * Float32Array.BYTES_PER_ELEMENT;
  const stateBuffer0 = device.createBuffer({
    size: stateBufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "state0",
  });
  const stateBuffer1 = device.createBuffer({
    size: stateBufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "state1"
  });
  const stateBuffer2 = device.createBuffer({
    size: stateBufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "state2"
  });
  // wave speeds
  const cBuffer = device.createBuffer({
    size: stateBufferSize,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "waveSpeed"
  });
  // cross section brightness buffer
  const brightnessBuffer = device.createBuffer({
    size: height * Float32Array.BYTES_PER_ELEMENT,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "crossSectionBrightness"
  });
  // time for wave generator
  const timeBuffer = device.createBuffer({
    size: Float32Array.BYTES_PER_ELEMENT,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST,
    label: "timeBuffer"
  });

  // wave speed array for updating
  const cArray = new Float32Array(numCells).fill(initC);
  // for undo
  const cArrPrev = new Float32Array(numCells).fill(initC);

  // Uniform Buffer
  // Layout: [width, height, dt, display, time, boundaryAbsorb, wavelength, amp, type(0=plane,1=point), waveOn, brightnessCrossSection, brightnessFilterStrength]
  const Uwidth = 0;
  const Uheight = 1;
  const Udt = 2;
  const Udisplay = 3;
  // const Utime = 4;
  const UboundaryAbsorb = 5;
  const Uwavelength = 6;
  const Uamp = 7;
  const Utype = 8;
  const UwaveOn = 9;
  const UbrightnessCrossSection = 10;
  const UbrightnessFilterStrength = 11;
  const uniformData = new Float32Array([width, height, dt,
    displayType[ui.displayType.value],
    0,
    ui.boundaryAbsorb.checked ? 1 : 0,
    parseInt(ui.wavelengthSlider.value),
    parseFloat(ui.amplitudeSlider.value),
    ui.planeWave.checked ? 0 : 1,
    1,
    width - 2,
    200
  ]);
  const uniformBuffer = device.createBuffer({
    size: uniformData.byteLength,
    usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    label: "uniform"
  });
  device.queue.writeBuffer(uniformBuffer, 0, uniformData.buffer);

  // Simulation Initialization

  function writeWaveState(f32offset, array) {
    const byteOffset = f32offset * Float32Array.BYTES_PER_ELEMENT;
    device.queue.writeBuffer(stateBuffer0, byteOffset, array.buffer);
    device.queue.writeBuffer(stateBuffer1, byteOffset, array.buffer);
    device.queue.writeBuffer(stateBuffer2, byteOffset, array.buffer);
  }

  function resetState() {
    const initialState = new Float32Array(numCells).fill(initU);

    writeWaveState(0, initialState);

    // reset time
    device.queue.writeBuffer(timeBuffer, 0, new Float32Array([0]));
  }

  function resetCBuffer() {
    cArray.fill(initC);
    device.queue.writeBuffer(cBuffer, 0, cArray.buffer);
  }

  function initializeState() {
    resetState();

    cArray.fill(initC);

    slitInterference();

    device.queue.writeBuffer(cBuffer, 0, cArray.buffer);

    device.queue.writeBuffer(brightnessBuffer, 0, (new Float32Array(height)).fill(0).buffer);
  }
  initializeState();

  const computeShaderCode = `
    struct Uniforms {
      width: f32,
      height: f32,
      dt: f32,
      display: f32,
      slot0: f32, // separate uniform for time
      boundaryAbsorb: f32,
      wavelength: f32,
      amp: f32,
      waveType: f32,
      waveOn: f32,
      crossSectionX: f32,
      brightnessFilterStrength: f32
    };
    @group(0) @binding(0) var<storage, read> statePrev: array<f32>;
    @group(0) @binding(1) var<storage, read> stateNow: array<f32>;
    @group(0) @binding(2) var<storage, read_write> stateNext: array<f32>;
    @group(0) @binding(3) var<storage, read> waveSpeed: array<f32>;
    @group(0) @binding(4) var<storage, read_write> brightness: array<f32>;
    @group(0) @binding(5) var<uniform> uniforms: Uniforms;
    @group(0) @binding(6) var<storage, read_write> time: f32;

    @compute @workgroup_size(16, 16)
    fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
      let pi2 = 6.28f;
      let x = i32(global_id.x);
      let y = i32(global_id.y);
      let width = i32(uniforms.width);
      let height = i32(uniforms.height);
      let index = y * width + x;
      let cdt = waveSpeed[index] * uniforms.dt;
      if (index == 0) {time += uniforms.dt;}

      if ((x >= width || y >= height) || (uniforms.boundaryAbsorb == 0 && (y >= height - 1 || y == 0)) || cdt <= 0) {
        return;
      }

      let frac = (cdt - 1) / (cdt + 1);

      let west = index - 1;
      let east = index + 1;
      let north = index - width;
      let south = index + width;
      let ne = north + 1;
      let nw = north - 1;
      let se = south + 1;
      let sw = south - 1;

      // 5 point stencil
      // let laplacian = (stateNow[east] + stateNow[west] + stateNow[south] + stateNow[north]
      // - 4 * stateNow[index]);

      // 9 point stencil
      let laplacian = (4 * (stateNow[east] + stateNow[west] + stateNow[south] + stateNow[north])
          + (stateNow[ne] + stateNow[nw] + stateNow[se] + stateNow[sw])
          - 20 * stateNow[index]) / 6.0f;

      stateNext[index] = 2 * stateNow[index] - statePrev[index] + cdt * cdt * laplacian;

      // absorbing boundaries - Mur's first order - right/left and c<0 barrier
      // x boundaries
      if (x == width - 1 || waveSpeed[east] < 0) {
        stateNext[index] = stateNow[west] + frac * (stateNext[west] - stateNow[index]);
        // return;
      }
      if (x == 0 || waveSpeed[west] < 0) {
        stateNext[index] = stateNow[east] + frac * (stateNext[east] - stateNow[index]);
        // return;
      }
      
      // y boundaries - only when enabled
      if ((uniforms.boundaryAbsorb == 1 && (y >= height - 1 || y <= 0)) || waveSpeed[north] < 0 || waveSpeed[south] < 0) {
        if (y == height - 1 || waveSpeed[south] < 0) {
          stateNext[index] = stateNow[north] + frac * (stateNext[north] - stateNow[index]);
          // second order mur abc
          // stateNext[index] = -statePrev[north] + frac * (stateNext[north] + statePrev[index])
          //   + (2.0 * (stateNow[index] + stateNow[north])
          //     + cdt * cdt / 2 * (stateNow[east] + stateNow[west] + stateNow[nw] + stateNow[ne] - 2 * (stateNow[index] + stateNow[north]))
          //   ) / (cdt + 1);
        }
        if (y == 0 || waveSpeed[north] < 0) {
          stateNext[index] = stateNow[south] + frac * (stateNext[south] - stateNow[index]);
        }
        return;
      }

      // Wave generators
      if (uniforms.waveOn == 1) {
        let waveGen = uniforms.amp * sin(time * pi2 / uniforms.wavelength);
        if (uniforms.waveType == 0 && x == 1) {
          // plane wave
          stateNext[index] = waveGen;
        } else if (uniforms.waveType == 1 && x == 50 && y == height / 2) {
          // point source
          stateNext[index] = 10 * waveGen;
        }
      }

      if (uniforms.display == 2 && x == i32(uniforms.crossSectionX)) {
        // brightness[y] = brightness[y] * uniforms.brightnessFilterStrength + stateNow[index] * stateNow[index] / 10;
        brightness[y] += (stateNow[index] * stateNow[index] - brightness[y]) / uniforms.brightnessFilterStrength; // low pass filter
        // brightness[y] += stateNow[index] * stateNow[index] / 10;
      } else if (uniforms.display < 2) {
        brightness[y] = 0;
      }
    }
  `;
  const computeModule = device.createShaderModule({ code: computeShaderCode, label: "computeModule" });
  const computePipeline = device.createComputePipeline({
    layout: 'auto',
    compute: { module: computeModule, entryPoint: 'main' },
    label: "computePipeline"
  });
  const computeBindGroup = (bufferPrev, bufferNow, bufferNext) => device.createBindGroup({
    layout: computePipeline.getBindGroupLayout(0),
    entries: [
      { binding: 0, resource: { buffer: bufferPrev } },
      { binding: 1, resource: { buffer: bufferNow } },
      { binding: 2, resource: { buffer: bufferNext } },
      { binding: 3, resource: { buffer: cBuffer } },
      { binding: 4, resource: { buffer: brightnessBuffer } },
      { binding: 5, resource: { buffer: uniformBuffer } },
      { binding: 6, resource: { buffer: timeBuffer } }
    ],
    label: "computeBindGroup"
  });

  // Render Pipeline (Visualization)
  const renderShaderCode = `
    struct Uniforms {
      width: f32,
      height: f32,
      dt: f32,
      display: f32,
      time: f32,
      boundaryAbsorb: f32,
      wavelength: f32,
      amp: f32,
      waveType: f32,
      waveOn: f32,
      crossSectionX: f32,
      brightnessFilterStrength: f32
    };
    @group(0) @binding(0) var<storage, read> state: array<f32>;
    @group(0) @binding(1) var<uniform> uniforms: Uniforms;
    @group(0) @binding(2) var<storage, read> waveSpeed: array<f32>;
    @group(0) @binding(3) var<storage, read> brightness: array<f32>;

    struct VertexOut {
      @builtin(position) position: vec4<f32>,
      @location(0) fragCoord: vec2<f32>,
    };

    @vertex
    fn vs_main(@builtin(vertex_index) vertexIndex: u32) -> VertexOut {
      var pos = array<vec2<f32>, 3>(
        vec2<f32>(-1.0,  3.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0, -1.0),
      );
      var output: VertexOut;
      output.position = vec4<f32>(pos[vertexIndex], 0.0, 1.0);
      output.fragCoord = 0.5 * (pos[vertexIndex] + vec2<f32>(1.0)) * vec2<f32>(uniforms.width, uniforms.height);
      return output;
    }

    fn colorMap(value:f32) -> vec3<f32> {
      //return vec3<f32>(value, 1.0 - abs(value - 0.5), 1.0 - value); // rgb
      return vec3<f32>(value, abs(value) / 2 - 0.5, -value); // yrkbc
    }

    @fragment
    fn fs_main(@location(0) fragCoord: vec2<f32>) -> @location(0) vec4<f32> {
      let x = i32(fragCoord.x);
      let y = i32(fragCoord.y);
      let width = i32(uniforms.width);
      let index = y * width + x;
      
      let west = index - 1;
      let east = index + 1;
      let north = index - width;
      let south = index + width;
      
      // Render barriers and large changes in c
      if ((waveSpeed[index] <= 0
           || abs(waveSpeed[north] - waveSpeed[south]) > 0.02
           || abs(waveSpeed[east] - waveSpeed[west]) > 0.02
          ) && uniforms.display <= 1
        ) {
        return vec4<f32>(1.0);
      }
      // Display wave state
      if (uniforms.display == 0) {
        // normal wave display
        return vec4<f32>(colorMap(state[index]), 1.0);
      } else if (uniforms.display == 1) {
        // wave speed display
        return vec4<f32>(colorMap(waveSpeed[index] - 1), 1.0);
      } else {
        // end screen display
        return vec4<f32>(brightness[y], brightness[y], brightness[y], 1.0);
        // return vec4<f32>(vec3<f32>(brightness[y] / uniforms.time), 1.0);
        // return vec4<f32>(colorMap(state[(y + 1) * width - 1]), 1.0);
      }
    }
  `;
  const renderModule = device.createShaderModule({ code: renderShaderCode, label: "renderModule" });
  const renderPipeline = device.createRenderPipeline({
    layout: 'auto',
    vertex: { module: renderModule, entryPoint: 'vs_main' },
    fragment: { module: renderModule, entryPoint: 'fs_main', targets: [{ format: swapChainFormat }] },
    primitive: { topology: 'triangle-list' },
    label: "renderPipeline"
  });
  const renderBindGroup = (stateBufferForRender) => device.createBindGroup({
    layout: renderPipeline.getBindGroupLayout(0),
    entries: [
      { binding: 0, resource: { buffer: stateBufferForRender } },
      { binding: 1, resource: { buffer: uniformBuffer } },
      { binding: 2, resource: { buffer: cBuffer } },
      { binding: 3, resource: { buffer: brightnessBuffer } },
    ],
    label: "renderBindGroup"
  });

  // Presets


  /**
   * Updates the wave speed of a cell
   * @param x x coordinate
   * @param y y coordinate
   * @param value new value of c
   */
  function updateCell(x, y, value) {
    if (x >= 0 && x < width && y >= 0 && y < height) {
      const index = y * width + x;
      cArray[index] = value;
      if (value <= 0) {
        writeWaveState(index, new Float32Array([initU]));
      }
    }
  }

  /**
   * Slit interference
   * @param slitWidth width of the slit
   * @param spacing spacing between slit centers
   * @param offsetLeft distance from left boundary to center of lens
   */
  function slitInterference(slitWidth = 20, spacing = 200, nSlits = 2, offsetLeft = 200) {
    const halfSpacing = Math.round(spacing / 2);
    const radius = Math.round(slitWidth / 2);
    
    for (let i = 0; i < halfHeight; i++) {
      // symmetric around half height
      if ((Math.abs((i + ((nSlits % 2 == 0) ? 0 : halfSpacing)) % spacing - halfSpacing) > radius) // fill outside of holes
        || (i > radius + (halfSpacing - 1) * nSlits) // fill the rest of the barrier
      ) {
        updateCell(offsetLeft, halfHeight + i, -1);
        updateCell(offsetLeft, halfHeight - i, -1);
      }
    }
  }


  /**
   * Lens generator
   * @param type Type configuration : "elliptical/parabolic" + "normal/fresnel"
   * @param r height of the lens
   * @param aspectRatio ratio of height to width
   * @param dir direction of curve : -1 left, 0 symmetrical, 1 right
   * @param convex convex (true) or concave (false)
   * @param depthRatio ratio of fresnel ridge height to half thickness, for fresnel only
   * @param refractiveIndex ratio of external wave speed to internal wave speed
   * @param offsetLeft distance from left simulation boundary to center of lens
   */
  function lensGenerator(type = "ellipticalnormal", r = height / 4, aspectRatio = 4, dir = 0, convex = true, depthRatio = 5, refractiveIndex = 1.2, offsetLeft = 300) {
    r = Math.ceil(r);
    dir *= convex ? 1 : -1;
    const yLimit = type.includes("fresnel") ? halfHeight : r;
    const lensWidth = Math.ceil(r / aspectRatio);
    const newC = initC / refractiveIndex;
    const a = aspectRatio * r;
    const depthRatioInv = 1 / depthRatio;
    const fresnelOffset = convex ? (10 - lensWidth * (1 - depthRatioInv)) : 0;

    // equations for lens shapes
    const lensShapes = Object.freeze({
      ellipticalnormal: (i, j) => (i * i + (aspectRatio * j) ** 2) < r * r,
      ellipticalfresnel: (i, j) => (i * i) % (r * r * (1 - (1 - depthRatioInv) ** 2)) + (aspectRatio * (j - fresnelOffset)) ** 2 < r * r,
      parabolicnormal: (i, j) => (i * i / a) + j - (convex ? 20 : 0) < lensWidth,
      parabolicfresnel: (i, j) => (i * i / a) % (lensWidth * depthRatioInv) + (j - fresnelOffset) < lensWidth,
    });

    for (let i = -yLimit; i < yLimit; i++) {
      for (let j = 0; j < lensWidth + 20; j++) {
        const offset = (convex ? 0 : lensWidth + 10);
        if (convex === lensShapes[type](i, j)) {
          if (dir >= 0) updateCell(j + offsetLeft - offset, i + halfHeight, newC);
          if (dir <= 0) updateCell(-j + offsetLeft + offset, i + halfHeight, newC);
        }
      }
    }
  }

  /**
   * Parabolic reflector generator
   * @param a width factor
   * @param dir direction of concavity: 1 = concave, -1 = convex to incoming waves
   * @param offsetLeft distance from left boundary to center of lens
   */
  function parabolicReflector(a = 1000, dir = 1, offsetLeft = width - 400) {
    a = Math.round(a / 25);
    for (let i = -halfHeight; i < halfHeight; i++) {
      for (let j = -halfWidth; j < halfWidth; j++) {
        if (Math.abs((i / a) ** 2 + dir * j) < 32 / a)
          updateCell(j + offsetLeft, i + halfHeight, 0);
      }
    }
  }
  /**
   * Circular reflector generator
   * @param r radius
   * @param dir direction of concavity: 1 = concave, -1 = convex to incoming waves
   * @param offsetLeft distance from left boundary to center of lens
   */
  function circularReflector(r = halfHeight, dir = 1, offsetLeft = Math.ceil(3 * width / 4)) {
    r = Math.ceil(r);
    for (let i = -r; i <= r; i++) {
      for (let j = -r; j <= r; j++) {
        if (Math.abs(i ** 2 + j ** 2 - r ** 2) <= r && dir * j >= 0) {
          updateCell(j - dir * r + offsetLeft, i + halfHeight, 0);
        }
      }
    }
  }


  // Settings and interaction
  {
    ui.cClear.addEventListener("click", () => {
      resetCBuffer();
      activePreset = null;
    });
    ui.reInit.addEventListener("click", () => resetState());

    ui.displayType.addEventListener("change", () => {
      const mode = ui.displayType.value;
      const modeVal = mode === "wave" ? 0 : mode === "speed" ? 1 : 2;
      uniformData[Udisplay] = modeVal;
      // device.queue.writeBuffer(uniformBuffer, 0, uniformData.buffer);
    });

    // dt
    ui.dt.addEventListener("input", () => {
      ui.dtValue.textContent = dt = uniformData[Udt] = parseFloat(ui.dt.value);
    });

    // Simulation speed
    ui.simSpeedSlider.addEventListener("input", () => {
      ui.simSpeedValue.textContent = dtPerFrame = parseInt(ui.simSpeedSlider.value);
    });

    // ABC condition toggle
    ui.boundaryAbsorb.addEventListener("click", () => {
      uniformData[UboundaryAbsorb] = ui.boundaryAbsorb.checked ? 1 : 0;
    });

    // Wave generator wavelength
    ui.wavelengthSlider.addEventListener("input", () => {
      uniformData[Uwavelength] = ui.wavelengthValue.textContent = parseInt(ui.wavelengthSlider.value);
    });
    // Wave generator amplitude
    ui.amplitudeSlider.addEventListener("input", () => {
      uniformData[Uamp] = ui.amplitudeValue.textContent = ampVal = parseFloat(ui.amplitudeSlider.value);
    });

    // Wave generator types
    ui.planeWave.addEventListener("input", () => {
      uniformData[Utype] = 0;
      resetState();
    });
    ui.pointSource.addEventListener("input", () => {
      uniformData[Utype] = 1;
      resetState();
    });

    // Toggle wave generator
    ui.waveToggle.addEventListener("click", () => {
      waveOn = !waveOn;

      if (waveOn) {
        device.queue.writeBuffer(timeBuffer, 0, new Float32Array([0]));
        uniformData[UwaveOn] = 1;
        if (uniformData[Uamp] < ampVal) {
          uniformData[Uamp] = ampVal;
        }
      }
    });

    // Preset settings
    function updatePreset() {
      resetCBuffer(); // disable to add
      const ol = parseInt(ui.offsetLeftSlider.value);
      switch (activePreset) {
        case Preset.slits:
          slitInterference(
            parseInt(ui.dsWidthSlider.value),
            parseInt(ui.dsSpacingSlider.value),
            parseInt(ui.dsNumSlits.value),
            ol
          );
          break;
        case Preset.lens:
          const data = new FormData(ui.lensOptions);
          const options = {};
          for (entry of data) {
            options[entry[0]] = entry[1];
          }
          const config = options["lensShape"] + options["lensType"];
          const r = parseInt(ui.lensRadiusSlider.value);
          const ar = parseFloat(ui.lensARSlider.value);
          const dir = parseInt(ui.lensDir.value);
          const curve = options.lensCurve == "convex";
          const ir = parseFloat(ui.lensIRSlider.value);
          const dr = parseFloat(ui.lensDRSlider.value);

          lensGenerator(config, r, ar, dir, curve, dr, ir, ol);
          break;
        case Preset.circularReflector:
          circularReflector(
            parseInt(ui.rRadiusSlider.value),
            parseInt(ui.rDir.checked ? -1 : 1),
            ol
          );
          break;
        case Preset.parabolicReflector:
          parabolicReflector(
            parseInt(ui.rRadiusSlider.value),
            parseInt(ui.rDir.checked ? -1 : 1),
            ol
          );
      }
      device.queue.writeBuffer(cBuffer, 0, cArray);
    }

    ui.preset.addEventListener("input", () => {
      ui.dsSettings.classList.add("hidden");
      ui.lensSettings.classList.add("hidden");
      ui.rSettings.classList.add("hidden");
      switch (Preset[ui.preset.value]) {
        case Preset.slits:
          ui.dsSettings.classList.remove("hidden");
          break;
        case Preset.lens:
          ui.lensSettings.classList.remove("hidden");
          break;
        case Preset.circularReflector:
          ui.rSettings.classList.remove("hidden");
          break;
        case Preset.parabolicReflector:
          ui.rSettings.classList.remove("hidden");
          break;
      }
    });

    HTMLElement.prototype.addOutput = function (valueOut, fixedPrecision = 0) {
      this.addEventListener("input", () => {
        if (valueOut) valueOut.textContent = fixedPrecision > 0 ? parseFloat(this.value).toFixed(fixedPrecision) : parseInt(this.value);
        updatePreset();
      });
    }

    ui.offsetLeftSlider.addOutput(ui.offsetLeftValue);
    ui.dsSpacingSlider.addOutput(ui.dsSpacingValue);
    ui.dsWidthSlider.addOutput(ui.dsWidthValue);
    ui.dsNumSlits.addOutput();
    ui.lensRadiusSlider.addOutput(ui.lensRadiusValue);
    ui.lensARSlider.addOutput(ui.lensARValue, 1);
    ui.lensIRSlider.addOutput(ui.lensIRValue, 2);
    ui.lensDir.addOutput();
    ui.lensDRSlider.addOutput(ui.lensDRValue, 1);
    ui.rRadiusSlider.addOutput(ui.rRadiusValue);

    ui.loadPreset.addEventListener("click", () => {
      activePreset = Preset[ui.preset.value];
      if (activePreset === Preset.circularReflector || activePreset === Preset.parabolicReflector)
        ui.offsetLeftSlider.value = ui.offsetLeftValue.textContent = halfWidth;
      updatePreset();
    });
  }


  // Barrier Setting (Mouse Click)
  {
    let isDrawing = false;
    let newC = null;
    canvas.addEventListener("mousedown", (event) => {
      if (uniformData[Udisplay] < 2) {
        isDrawing = true;
        placeBarrier(event);
      }
    });
    canvas.addEventListener("pointermove", (event) => {
      if (isDrawing) {
        event.getCoalescedEvents().forEach((event) => placeBarrier(event));
      }
    });
    canvas.addEventListener("mouseup", () => {
      isDrawing = false;
      newC = null;
    });
    function placeBarrier(event) {
      const rect = canvas.getBoundingClientRect();
      const x = Math.floor(event.clientX - rect.left);
      const y = height - Math.floor(event.clientY - rect.top);
      const index = y * width + x;
      if (newC === null) newC = cArray[index] === initC ? (ui.drawnBarrierAbsorb.checked ? -1 : 0) : initC;
      cArray[index] = newC;
      const indexOffsets = [0, 1, -1, width, -width, width + 1, width - 1, -width + 1, -width - 1];
      indexOffsets.forEach((offset) => {
        const totalOffset = index + offset;
        device.queue.writeBuffer(cBuffer, totalOffset * Float32Array.BYTES_PER_ELEMENT, new Float32Array([newC]))
        writeWaveState(totalOffset, new Float32Array([initU]));
      });
    }

    // Image Upload & Processing
    ui.blackRefractSlider.addOutput(ui.blackRefractValue, 2);
    ui.whiteRefractSlider.addOutput(ui.whiteRefractValue, 2);
    ui.imageApply.addEventListener("click", () => {
      if (!ui.imageUpload.files || ui.imageUpload.files.length === 0) return;
      resetCBuffer();
      const file = ui.imageUpload.files[0];
      const reader = new FileReader();
      reader.onload = (e) => {
        const img = new Image();
        img.onload = () => {
          // Get the UI scale factor.
          const uiScale = parseFloat(ui.imageScale.value);
          // Compute the maximum scale factor to fit the canvas.
          const fitScale = Math.min(width / img.width, height / img.height);
          // Final target size: original image scaled by fitScale and then by uiScale.
          const targetWidth = Math.round(img.width * fitScale * uiScale);
          const targetHeight = Math.round(img.height * fitScale * uiScale);
          // Create offscreen canvas to draw the scaled image.
          const offCanvas = document.createElement("canvas");
          offCanvas.width = targetWidth;
          offCanvas.height = targetHeight;
          const offCtx = offCanvas.getContext("2d");
          offCtx.drawImage(img, 0, 0, targetWidth, targetHeight);
          const imageData = offCtx.getImageData(0, 0, targetWidth, targetHeight);
          // Compute offsets to center the barrier image in the simulation grid.
          const offsetX = Math.floor((width - targetWidth) / 2);
          const offsetY = Math.floor((height - targetHeight) / 2);
          for (let j = 1; j < targetHeight; j++) { // unstable with j = 0
            for (let i = 0; i < targetWidth; i++) {
              const idx = ((targetHeight - j) * targetWidth + i) * 4;
              // Compute normalized brightness (average of R, G, B).
              const b = (imageData.data[idx] + imageData.data[idx + 1] + imageData.data[idx + 2]) / (3 * 255);
              const brightness = ui.barrierInvert.checked ? 1 - b : b;
              const simX = offsetX + i;
              const simY = offsetY + j;
              if (simX >= 0 && simX < width && simY >= 0 && simY < height) {
                updateCell(
                  simX,
                  simY,
                  (ui.blockBottom.checked && j === 1 || ui.blockTop.checked && j === targetHeight - 1)
                    ? -1
                    : (1 - brightness) * initC / parseFloat(ui.blackRefractSlider.value)
                    + brightness * initC / parseFloat(ui.whiteRefractSlider.value)
                );
              }
            }
          }
          device.queue.writeBuffer(cBuffer, 0, cArray.buffer);
        };
        img.src = e.target.result;
      };
      reader.readAsDataURL(file);
    });
  }


  // Simulation Loop

  const filterStrength = 10;
  let rafId;
  let jsTime = 0, lastFrameTime = performance.now(), deltaTime = 10, fps = 0;

  let [bufferPrev, bufferNow, bufferNext] = [stateBuffer0, stateBuffer1, stateBuffer2];

  function frame() {
    // now *= 1e-3;
    const startTime = performance.now();
    deltaTime += (startTime - lastFrameTime - deltaTime) / filterStrength;
    fps += (1e3 / deltaTime - fps) / filterStrength;
    lastFrameTime = startTime;


    // ease the wave off when disabled
    if (!waveOn && uniformData[Uamp] > 0) {
      uniformData[Uamp] -= decayRate * ampVal; //Math.min(0.2, decayRate / uniformData[Uwavelength] * ampVal * dtPerFrame);
      if (uniformData[Uamp] <= 0) {
        uniformData[Uamp] = 0;
        uniformData[UwaveOn] = 0;
      }
    }

    device.queue.writeBuffer(uniformBuffer, 0, uniformData.buffer);

    const commandEncoder = device.createCommandEncoder();
    // run several timesteps per frame
    for (let i = 0; i < dtPerFrame; i++) {
      const computePass = commandEncoder.beginComputePass();
      computePass.setPipeline(computePipeline);
      computePass.setBindGroup(0, computeBindGroup(bufferPrev, bufferNow, bufferNext));
      computePass.dispatchWorkgroups(Math.ceil(width / 16), Math.ceil(height / 16));
      computePass.end();

      // rotate the buffers by 1 time step
      [bufferPrev, bufferNow, bufferNext] = [bufferNow, bufferNext, bufferPrev];
    }
    const textureView = context.getCurrentTexture().createView();
    const renderPass = commandEncoder.beginRenderPass({
      colorAttachments: [{
        view: textureView,
        clearValue: { r: 0, g: 0, b: 0, a: 1 },
        loadOp: 'clear',
        storeOp: 'store',
      }],
    });
    renderPass.setPipeline(renderPipeline);
    renderPass.setBindGroup(0, renderBindGroup(bufferNow));
    renderPass.draw(3, 1, 0, 0);
    renderPass.end();

    device.queue.submit([commandEncoder.finish()]);

    jsTime += (performance.now() - startTime - jsTime) / filterStrength;

    //setTimeout(frame, 500);
    rafId = requestAnimationFrame(frame);
  }
  rafId = requestAnimationFrame(frame);

  const intID = setInterval(() => {
    ui.fps.textContent = fps.toFixed(1);
    ui.jsTime.textContent = jsTime.toFixed(2);
    ui.frameTime.textContent = deltaTime.toFixed(1);
    // ui.gpuTime.textContent = deltaTime - jsTime;
  }, 100);


  // Other event listeners
  ui.collapse.onclick = () => {
    ui.collapse.innerText = ui.collapse.innerText === ">" ? "<" : ">";
    if (ui.panel.classList.contains("hidden")) {
      ui.panel.classList.remove("hidden");
      ui.collapse.classList.remove("inactive");
    } else {
      ui.panel.classList.add("hidden");
      ui.collapse.classList.add("inactive");
    }
  };
  
  window.onresize = () => {
    exit = true;
    clearInterval(intID);
    cancelAnimationFrame(rafId);
    main();
  };
}

main();