////////////////////////////////////////////////////////////////////////////////
//
// xrscene.js
// Copyright (c) 2021 Matthew Arcus
// MIT License (MIT)
//
// Set up WebXR fragment shader, either in VR or AR mode.
//
////////////////////////////////////////////////////////////////////////////////

"use strict"

const VR_MODE = 0;
const AR_MODE = 1;

function xrscene(mode,xrButton) {
    if (mode != VR_MODE && mode != AR_MODE) {
        throw new Error('Invalid XR mode')
    }
    let sessiontype = null;
    let modename = null;
    if (mode == VR_MODE) {
        sessiontype = 'immersive-vr';
        modename = 'VR';
    } 
    if (mode == AR_MODE) {
        sessiontype = 'immersive-ar';
        modename = 'AR';
    } 

    // Experiment with 'local-floor', 'unbounded', etc.
    // 'local' works on my phone
    let referencespace = 'local';

    // If shaderfile is defined, load locally,
    let shaderfile = null
    //shaderfile = "goursat.glsl";
    shaderfile = "pentagram.glsl";

    // else try and load a Shadertoy shader
    let shaderID = "3d2GDt" // Shadertoy Goursat
    //let shaderID = "4sX3Rn" // iq's menger sponge
    //let shaderID = "XdGczw"     // parallepiped
    //let shaderID = "4tSBDz"     // inverted spheres

    let xrSession = null;
    let xrRefSpace = null;
    let gl = null;
    let renderer = null;
    let framecount = 0;
    let selectcount = 0;

    const devicePixelRatio = window.devicePixelRatio || 1; // This is 3 on my phone
    //if (devicePixelRatio != 1) alert("devicePixelRatio: " + devicePixelRatio); // Find out DPR

    function makeShader(source, shadertype) {
        const shader = gl.createShader(shadertype);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        const infolog = gl.getShaderInfoLog(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            if (infolog) alert("Info for " + shadertype + " shader: " + infolog);
            return shader;
        } else {
            let typestring = "<unknown>";
            if (shadertype == gl.VERTEX_SHADER) typestring = "vertex";
            else if (shadertype == gl.FRAGMENT_SHADER) typestring = "fragment";
            alert("Error for " + shadertype + " shader: " + infolog);
        }
    }

    function Renderer(fsource) {
        this.transformMatrix = mat4.create(); // Cache the transform matrix here

        // Compile and link the shaders.
        const vsource = [
            "#version 300 es",
            "in vec3 aVertexPosition;",
            "out vec2 vTextureCoord;",
            "void main(void) {",
            "  gl_Position = vec4(aVertexPosition,1.0);",
            "  vTextureCoord = gl_Position.xy;",
            "}"
        ].join("\n");

        const preamble = [
            "#version 300 es",
            "precision highp float;",
            // XR uniforms
            // Declare projection & view in shader, if needed.
            //"uniform mat4 iProjection;", 
            //"uniform mat4 iView;",
            "uniform mat4 transformMatrix;",

            // A partial set of Shadertoy uniforms
            "uniform float iTime;",
            "uniform int iFrame;",
            "uniform vec4 iResolution;",
            "uniform vec4 iMouse;",
            "out vec4 outColor;",
            "in vec2 vTextureCoord;",

            "void mainVR(out vec4 color, vec2 fragCoord, vec3 p, vec3 q);",
            "void main() {",
            "  vec4 screenpos = vec4(vTextureCoord,0,1);", // The "screen position", -1 <= z <= 1
            "  vec4 eye = vec4(0,0,1,0);",         // z-infinity
            "  eye = transformMatrix*eye;",
            "  screenpos = transformMatrix*screenpos;",
            "  vec3 p = eye.xyz/eye.w;",
            "  vec3 r = screenpos.xyz/screenpos.w;",
            "  mainVR(outColor,vTextureCoord,p,normalize(r-p));",
            (mode == VR_MODE ? "  outColor.a = 1.0;" : ""), // Use alpha blending for AR
            "}",
            "#define HW_PERFORMANCE 0", // Some Shadertoys use this to set AA
            "#line 1",
            ""
        ].join("\n");

        if (!shaderfile) {
            // Pull data out from Shadertoy JSON response.
            // This is probably a bit fragile.
            const json = JSON.parse(fsource);
            let error = json['Error'];
            if (error) throw new Error(error);
            const shader = json['Shader']
            console.log(shader.info);
            fsource = "";
            for (const pass of shader.renderpass) {
                if (pass.type == 'common') fsource = pass.code+fsource;
                if (pass.type == 'image') fsource = fsource+pass.code;
            }
            if (!fsource) throw new Error('No shader source found');
        }
        fsource = preamble + fsource;
        const vshader = makeShader(vsource,gl.VERTEX_SHADER);
        const fshader = makeShader(fsource,gl.FRAGMENT_SHADER);
        if (!vshader || !fshader) throw new Error("compilation failed");
        const program = gl.createProgram();
        if (!program) throw new Error("program creation failed");
        gl.attachShader(program, vshader);
        gl.attachShader(program, fshader);
        gl.linkProgram(program);
        gl.validateProgram(program); // Check all well
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            throw(new Error("Unable to initialize program: " + gl.getProgramInfoLog(program)));
        }

        // Two triangles fill the screen
        const vertices = [
            1.0, 1.0, 0.0,
            -1.0, 1.0, 0.0,
            1.0,-1.0, 0.0,
            -1.0,-1.0, 0.0
        ];
        this.vertBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
        this.program = program;
    }

    Renderer.prototype.draw = function (projectionMatrix, viewMatrix, timestamp) {
        if (this.start == null) this.start = timestamp;
        const transformMatrix = this.transformMatrix;
        const program = this.program;
        gl.useProgram(program);
        // transform = inverse(projection*view);
        mat4.mul(transformMatrix,projectionMatrix,viewMatrix);
        mat4.invert(transformMatrix,transformMatrix);
        const index = gl.getAttribLocation(program, "aVertexPosition");
        if (index >= 0) {
            gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
            gl.enableVertexAttribArray(index);
            gl.vertexAttribPointer(index, 3, gl.FLOAT, false, 3*4, 0*4);
        }
        gl.uniform1i(gl.getUniformLocation(program, "iFrame"), framecount++);
        gl.uniform1f(gl.getUniformLocation(program, "iTime"), (timestamp-this.start)/1000);
        // Parameter 2 is `transpose` - set it to false!
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "iView"), false, viewMatrix);
        //gl.uniformMatrix4fv(gl.getUniformLocation(program, "iProjection"), false, projectionMatrix);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "transformMatrix"), false, transformMatrix);
        // drawingBufferWidth/Height should be the right thing
        const width = gl.canvas.drawingBufferWidth, height = gl.canvas.drawingBufferHeight;
        gl.uniform4f(gl.getUniformLocation(program, "iResolution"),
                     width*devicePixelRatio, height*devicePixelRatio,0,0);
        // Should pass in ray to centre of image instead of mouse.
        // "selectcount" is number of select events seen.
        gl.uniform4f(gl.getUniformLocation(program, "iMouse"),0,0,0,selectcount);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        const err = gl.getError();
        if (err) console.log("GL Error: ", err); // Alert? Exception?
    };

    // Called either when the user has explicitly ended the session by calling
    // session.end() or when the UA has ended the session for any reason.
    // At this point the session object is no longer usable and should be
    // discarded.
    function onSessionEnded(event) {
        console.log("WebXR session ended");
        xrSession = null;
        xrButton.textContent = 'Enter ' + modename

        // In this simple case discard the WebGL context too, since we're not
        // rendering anything else to the screen with it.
        gl = null;
    }
    
    function onXRFrame(timestamp, frame) {
        console.log("XRFrame");
        const session = frame.session;
        session.requestAnimationFrame(onXRFrame);
        const pose = frame.getViewerPose(xrRefSpace);
        if (false) {
            // Find input sources
            for (let inputSource of frame.session.inputSources) {
                let targetRayPose = frame.getPose(inputSource.targetRaySpace, xrRefSpace);
                if (!targetRayPose) {
                    continue;
                }
                alert(targetRayPose);
            }
        }
        if (pose) {
            const baseLayer = session.renderState.baseLayer;
            gl.bindFramebuffer(gl.FRAMEBUFFER, baseLayer.framebuffer);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            for (const view of pose.views) {
                const viewport = baseLayer.getViewport(view);
                gl.viewport(viewport.x, viewport.y, viewport.width, viewport.height);
                renderer.draw(view.projectionMatrix, view.transform.inverse.matrix, timestamp);
            }
        }
    }

    function onSelectStart(event) {
        //alert("onSelectStart");
    }
    function onSelectEnd(event) {
        //alert("onSelectEnd");
    }
    function onSelect(event) {
        //alert("onSelect");
        selectcount++;
    }
    
    function makeshaderurl(shaderID) {
        const qparam = new Date().getTime();  // Skip caching
        if (shaderfile) return shaderfile + "?" + qparam;
        else if (shaderID) return "https://www.shadertoy.com/api/v1/shaders/" + shaderID + "?key=fdntwh&" + qparam;
        else throw new Error('No shader specified');
    }

    function onSessionStarted (session) {
        xrSession = session;
        xrButton.textContent = 'Exit ' + modename;
        session.addEventListener('end', onSessionEnded);

        if (true) {
            session.addEventListener('selectstart', onSelectStart);
            session.addEventListener('selectend', onSelectEnd);
            // By listening for the 'select' event we can find out when the user has
            // performed some sort of primary input action and respond to it.
            session.addEventListener('select', onSelect);
        }

        const webglCanvas = document.createElement('canvas');
        gl = webglCanvas.getContext('webgl2', { xrCompatible: true })
        if (!gl) throw new Error('getContext failed');

        // Send off a request for the fragment shader code.
        const request = new XMLHttpRequest();
        request.open("GET", makeshaderurl(shaderID));
        request.onreadystatechange = function() {
            if (request.readyState === 4) {
                function setRefSpace(refSpace) {
                    xrRefSpace = refSpace;
                    session.requestAnimationFrame(onXRFrame);
                }
                renderer = new Renderer(request.responseText);
                session.updateRenderState({ baseLayer: new XRWebGLLayer(session, gl) });
                session.requestReferenceSpace(referencespace).
                    then(setRefSpace,
                         error => {
                             // viewer reference space should be defined
                             alert("defaulting to 'viewer' reference space");
                             session.requestReferenceSpace('viewer').then(setRefSpace);
                         });
            }
        }
        request.send(null); // No body
    }

    // Called when the user clicks the button to enter XR. If we don't have a
    // session we'll request one, and if we do have a session we'll end it.
    function onButtonClicked() {
        if (xrSession) xrSession.end();
        else navigator.xr.requestSession(sessiontype).then(onSessionStarted);
    }

    // Find browser capabilities
    if (navigator.xr) {
        navigator.xr.isSessionSupported(sessiontype).then((supported) => {
            if (supported) {
                xrButton.addEventListener('click', onButtonClicked);
                xrButton.textContent = 'Enter ' + modename;
                xrButton.disabled = false;
            }
        });
    }
}
