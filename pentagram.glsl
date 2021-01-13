#if defined HEADER
precision mediump float;
uniform float iTime;
uniform vec2 iResolution;
uniform vec4 iMouse;

void mainImage(out vec4 outColor, vec2 fragCoord);

out vec4 outColor;

void main() {
  mainImage(outColor, gl_FragCoord.xy);
  outColor.a = 1.0;
  //outColor = vec4(hsv(gl_FragCoord.x/iResolution.x),1.0);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// (c) Matthew Arcus 2017
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Glowing festive pentagram: ray trace 5 rotations of basic line, use
// distance of ray to line to determine color.
////////////////////////////////////////////////////////////////////////////////

float N = 5.0; // Number of lines
const float A = 0.6; // Light amplitude
const float K = 0.4; // Concentration
const float R = 1.0; // Radius
const float PI = 3.141592654;

float det(vec2 p, vec2 q) {
  return determinant(mat2(p,q));
}

vec2 closest(vec3 p,vec3 q,vec3 r,vec3 s) {
  // Use Cramer's rule to solve linear system
  // Assumes q and s are unit vectors
  // No cross products, 4 dot products, 3 2x2 determinants
  vec2 c0 = vec2(1.0,dot(q,s));
  vec2 c1 = vec2(-dot(q,s),-1.0);
  vec2 a = vec2(dot(r-p,q),dot(r-p,s));
  return vec2(det(a,c1),det(c0,a))/det(c0,c1);
}

vec2 closest1(vec3 p,vec3 q,vec3 r,vec3 s) {
  // Matrix inverse solution
  // Assumes q and s are unit vectors
  float k = dot(q,s);
  mat2 m = mat2(1,k,-k,1);
  return inverse(m)*vec2(dot(r-p,q),dot(r-p,s));
}

// Smooth HSV to RGB conversion 
// Function by iq, from https://www.shadertoy.com/view/MsS3Wc
vec3 h2rgb(float h) {
  vec3 rgb = clamp( abs(mod(h*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
  rgb = rgb*rgb*(3.0-2.0*rgb); // cubic smoothing	
  return rgb;
}

// Quaternion to rotation matrix
// Assumes normalized
mat3 qrot(vec4 q) {
  float x = q.x, y = q.y, z = q.z, w = q.w;
  float x2 = x*x, y2 = y*y, z2 = z*z;
  float xy = x*y, xz = x*z, xw = x*w;
  float yz = y*z, yw = y*w, zw = z*w;
  return 2.0*mat3(0.5-y2-z2, xy+zw, xz-yw,
                  xy-zw, 0.5-x2-z2, yz+xw,
                  xz+yw, yz-xw, 0.5-x2-y2);
}

vec2 rotate(vec2 p, float t) {
  return p * cos(t) + vec2(p.y, -p.x) * sin(t);
}

vec4 mainFun(vec3 p, vec3 q, vec3 rcentre) {
  // r+js is polygon line, to be rotated in loop
  // Rotation axis
  vec3 axis = normalize(vec3(1,1,cos(0.1*iTime)));
  float phi = iTime*0.15;
  mat3 n = qrot(vec4(sin(phi)*axis,cos(phi)));
  p = n*p; q = n*q;
  float scale = 2.0;
  p *= scale;
  q = normalize(q);
  float mindist = 1e10;
  vec4 color = vec4(0); // Accumulate color here
  vec3 r = vec3(0,1,0);
  N += floor(mod(iMouse.w,6.0));
  float theta = 2.0*PI/N;
  mat2 m = mat2(cos(theta),sin(theta),-sin(theta),cos(theta));
  for (float i = 0.0; i < N; i++, r.xy*=m) {
    r = normalize(r);
    vec3 s = vec3(r.y,-r.x,0);
    vec2 k = closest(p,q,r,s);
    if (k.x < 0.0) continue;
    float d = distance(p+k.x*q,r+k.y*s);
    float h = fract(0.3*(-iTime+log(1.0+abs(k.y))));
    vec4 basecolor = vec4(h2rgb(h),1.0);
    color += A*(1.0-pow(smoothstep(0.0,R,d),K))*basecolor;
  }
  color = sqrt(color);
  return color;
}

vec4 mainXR(vec3 eye, vec3 ray, vec3 rcentre) {
  // Should have a proper model transform here
  vec3 modeloffset = vec3(0,0,2.5);
  return mainFun(eye+modeloffset,ray,rcentre);
}

void mainVR( out vec4 fragColor, in vec2 fragCoord,
             in vec3 fragRayOrigin, in vec3 fragRayDir) {
  vec3 modeloffset = vec3(0,0,2.5);
  fragColor = mainFun(fragRayOrigin+modeloffset,fragRayDir,vec3(0));
}

vec3 transform(in vec3 p) {
  if (iMouse.x > 0.0) {
    float theta = (2.0*iMouse.y-iResolution.y)/iResolution.y*PI;
    float phi = (2.0*iMouse.x-iResolution.x)/iResolution.x*PI;
    p.yz = rotate(p.yz,theta);
    p.zx = rotate(p.zx,-phi);
  }
  return p;
}

void mainImage(out vec4 outColor, vec2 fragCoord) {
  vec2 xy = (2.0*fragCoord - iResolution.xy)/iResolution.y;
  // p+kq is viewing ray
  vec3 p = vec3(0,0,4);
  vec3 q = vec3(xy,-2);
  p = transform(p);
  q = transform(q);
  q = normalize(q);
  outColor = mainFun(p,q,vec3(0,0,-1));
}
