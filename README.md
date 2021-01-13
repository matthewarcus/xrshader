<link href='css/common.css' rel='stylesheet'></link>
<script src='js/gl-matrix-min.js'></script>
<script src='xrscene.js'></script>

# xrshader
Experimenting with fragment shaders for WebXR

WIP, but feedback welcome.

[Try It!](https://matthewarcus.github.io/xrshader/xrscene.html)

Tested with Chrome on Samsung A5 phone.
<button id="vr-button" class="barebones-button" disabled>VR not found</button>
<button id="ar-button" class="barebones-button" disabled>AR not found</button>
<script>xrscene(VR_MODE,document.getElementById('vr-button'))</script>
<script>xrscene(AR_MODE,document.getElementById('ar-button'))</script>
