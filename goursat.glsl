////////////////////////////////////////////////////////////////////////////////
//
// Goursat Quartic Surfaces
// Copyright (c) 2021 Matthew Arcus
// MIT License (MIT)
//
// http://mathworld.wolfram.com/GoursatsSurface.html
// https://www.mathcurve.com/surfaces/goursat/goursat.shtml
//
// Quartic surfaces with octahedral symmetry. Surface (including normals)
// is raytraced using analytic solution to quartic due to Lanczos and Kahan.
//
// XR version
//
////////////////////////////////////////////////////////////////////////////////

uniform mat4 iView; // XR view transform

// PARAMS
bool dorotate = true;

// Lighting
vec3 light = normalize(vec3(1,1,1));
float ambient = 0.4;
float diffuse = 0.6;
float specular = 0.8;
float specularpow = 10.0;
vec3 specularcolor = vec3(1);

// Debug
bool alert = false;

void assert(bool t) {
  if (!t) alert = true;
}

const float PI =  3.141592654;

vec2 rotate(vec2 p, float t) {
  return p * cos(t) + vec2(p.y, -p.x) * sin(t);
}

float sgn(float x) {
  return x < 0.0 ? -1.0: 1.0; // Return 1 for x == 0
}

float evalquadratic(float x, float A, float B, float C) {
  return (A*x+B)*x+C;
}

float evalcubic(float x, float A, float B, float C, float D) {
  return ((A*x+B)*x+C)*x+D;
}

// Quadratic solver from Kahan
int quadratic(float A, float B, float C, out vec2 res) {
  float b = -0.5*B, b2 = b*b;
  float q = b2 - A*C;
  if (q < 0.0) return 0;
  float r = b + sgn(b)*sqrt(q);
  if (r == 0.0) {
    res[0] = C/A;
    res[1] = -res[0];
  } else {
    res[0] = C/r;
    res[1] = r/A;
  }
  return 2;
}

// Numerical Recipes algorithm for solving cubic equation
int cubic(float a, float b, float c, float d, out vec3 res) {
  if (a == 0.0) {
    return quadratic(b,c,d,res.xy);
  }
  if (d == 0.0) {
    res.x = 0.0;
    return 1+quadratic(a,b,c,res.yz);
  }
  float tmp = a; a = b/tmp; b = c/tmp; c = d/tmp;
  // solve x^3 + ax^2 + bx + c = 0
  float Q = (a*a-3.0*b)/9.0;
  float R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  float R2 = R*R, Q3 = Q*Q*Q;
  if (R2 < Q3) {
    float X = clamp(R/sqrt(Q3),-1.0,1.0);
    float theta = acos(X);
    float S = sqrt(Q); // Q must be positive since 0 <= R2 < Q3
    res[0] = -2.0*S*cos(theta/3.0)-a/3.0;
    res[1] = -2.0*S*cos((theta+2.0*PI)/3.0)-a/3.0;
    res[2] = -2.0*S*cos((theta+4.0*PI)/3.0)-a/3.0;
    return 3;
  } else {
    float alpha = -sgn(R)*pow(abs(R)+sqrt(R2-Q3),0.3333);
    float beta = alpha == 0.0 ? 0.0 : Q/alpha;
    res[0] = alpha + beta - a/3.0;
    return 1;
  }
}

float qcubic(float B, float C, float D) {
  vec3 roots;
  int nroots = cubic(1.0,B,C,D,roots);
  // Sort into descending order
  if (nroots > 1 && roots.x < roots.y) roots.xy = roots.yx;
  if (nroots > 2) {
    if (roots.y < roots.z) roots.yz = roots.zy;
    if (roots.x < roots.y) roots.xy = roots.yx;
  }
  // And select the largest
  float psi = roots[0];
  // There _should_ be a positive root, but sometimes the cubic
  // solver doesn't find it directly (probably a double root
  // around zero).
  if (psi < 0.0) assert(evalcubic(psi,1.0,B,C,D) < 0.0);
  // If so, nudge in the right direction
  psi = max(1e-6,psi);
  // and give a quick polish with Newton-Raphson
  for (int i = 0; i < 3; i++) {
    float delta = evalcubic(psi,1.0,B,C,D)/evalquadratic(psi,3.0,2.0*B,C);
    psi -= delta;
  }
  return psi;
}

// The Lanczos quartic method
int lquartic(float c1, float c2, float c3, float c4, out vec4 res) {
  float alpha = 0.5*c1;
  float A = c2-alpha*alpha;
  float B = c3-alpha*A;
  float a,b,beta,psi;
  psi = qcubic(2.0*A-alpha*alpha, A*A+2.0*B*alpha-4.0*c4, -B*B);
  psi = max(0.0,psi);
  a = sqrt(psi);
  beta = 0.5*(A + psi);
  if (psi <= 0.0) {
    b = sqrt(max(beta*beta-c4,0.0));
  } else {
    b = 0.5*a*(alpha-B/psi);
  }
  int resn = quadratic(1.0,alpha+a,beta+b,res.xy);
  vec2 tmp;
  if (quadratic(1.0,alpha-a,beta-b,tmp) != 0) { 
    res.zw = res.xy;
    res.xy = tmp;
    resn += 2;
  }
  return resn;
}

int quartic(float A, float B, float C, float D, float E, out vec4 roots) {
  int nroots;
  // Solve for the smallest cubic term, this seems to give the least wild behaviour.
  if (abs(B/A) < abs(D/E)) {
    nroots = lquartic(B/A,C/A,D/A,E/A,roots);
  } else {
    nroots = lquartic(D/E,C/E,B/E,A/E,roots);
    for (int i = 0; i < nroots; i++) {
      roots[i] = 1.0/roots[i];
    }
  }
  assert(nroots == 0 || nroots == 2 || nroots == 4);
  return nroots;
}

struct Surface {
  vec4 params;
  vec3 p;
  int colorscheme;
};

// Equation: pp.pp + k(p.p)^2 + k'a^2(p.p) + k''a^4 = 0
// Derivative: 4ppp + 4k(p.p)p + 2k'a^2p
// Expansion with p => p+tr:
// pp => (p+tr)(p+tr) = pp + 2tpr + t^2rr
// pp.pp => (pp + 2tpr + t^2rr).(pp + 2tpr + t^2rr)
//  = pp.pp + 4tpp.pr + 6t^2pp.rr + 4t^3pr.rr + t^4rr.rr 
// p.p  => (p+tr).(p+tr) = p.p + 2tp.r + t^2r.r = p.p + 2tp.r + t^2
// (p.p)^2 = (p.p + 2tp.r + t^2)(p.p + 2tp.r + t^2) =
//         = p.p^2 + 4t^2(p.r)^2 + t^4 + 2(2t(p.p)(p.r) + (p.p)t^2 + 2t^3(p.r))
//         = p.p^2 + 4t^2(p.r)^2 + t^4 + 4t(p.p)(p.r) + 2(p.p)t^2 + 4t^3(p.r))
// ie.
// pp.pp + 4tpp.pr + 6t^2pp.rr + 4t^3pr.rr + t^4rr.rr +
// k(p.p^2 + t4[(p.p)(p.r)] + t^2[4(p.r)^2 + 2(p.p)] + t^3[4(p.r)] + t^4) +
// k'a^2(p.p + 2tp.r + t^2) +
// k''a^4
// collecting terms:
// t^0: pp.pp +   k(p.p)^2 +             k'a^2(p.p) + k''a^4
// t^1: 4pp.pr + 4k(p.p)(p.r) +         2k'a^2(p.r)
// t^2: 6pp.rr +  k[4(p.r)^2 + 2(p.p)] + k'a^2
// t^3: 4pr.rr + 4k(p.r)
// t^4: rr.rr +   k
//
// Can adjust p so that p.r = 0 & that simplifies things.

int goursatsurface(Surface surface, vec3 p, vec3 r, out vec4 roots) {
  float k = surface.params[0];
  float k1 = surface.params[1];
  float k2 = surface.params[2];
  float a = surface.params[3];
  vec3 pp = p*p;
  vec3 pr = p*r;
  vec3 rr = r*r;
  float p2 = dot(p,p);
  float a2 = a*a;
  float a4 = a2*a2;
#if 0
  // Unoptimized version
  float pdr = dot(p,r);
  float A = dot(rr,rr) + k;
  float B = 4.0*dot(pr,rr) + 4.0*k*pdr;
  float C = 6.0*dot(pp,rr) +     k*(4.0*pdr*pdr + 2.0*p2) + k1*a2;
  float D = 4.0*dot(pp,pr) + 4.0*k*(p2*pdr)           + 2.0*k1*a2*pdr;
  float E = dot(pp,pp)         + k*p2*p2 + k1*a2*p2 + k2*a4;
#else
  // Optimized, with p.r = 0
  float A =     dot(rr,rr) +     k;
  float B = 4.0*dot(pr,rr);
  float C = 6.0*dot(pp,rr) + 2.0*k*p2    + k1*a2;
  float D = 4.0*dot(pp,pr);
  float E =     dot(pp,pp) +     k*p2*p2 + k1*a2*p2 + k2*a4;
#endif
  return quartic(A,B,C,D,E,roots);
}

vec3 goursatnormal(Surface surface, vec3 p) {
  // 4ppp + 4k(p.p)p + 2k'a^2p
  float k = surface.params[0];
  float k1 = surface.params[1];
  float a = surface.params[3];
  return 4.0*p*p*p + 4.0*k*dot(p,p)*p + 2.0*k1*a*a*p;
}

vec3 applylighting(vec3 baseColor, vec3 p, vec3 n, vec3 r) {
  if (dot(r,n) > 0.0) n = -n; // Face forwards
  vec3 c = baseColor*ambient;
  c += baseColor*diffuse*(max(0.0,dot(light,n)));
  float s = pow(max(0.0,dot(reflect(light,n),r)),specularpow);
  c += specular*s*specularcolor;
  return c;
}

struct Result {
  vec3 p;
  vec3 n;
  vec3 basecolor;
  float t;
};

float gridline(vec3 p) {
  // Draw some gridlines on surface
  vec3 t = fract(p*4.0);
  t = min(t,1.0-t);
  float d = min(t.x,min(t.y,t.z));
  return smoothstep(0.02,0.025,d);
}

int dosurface(Surface surface, vec3 p0, vec3 r, out vec4 roots) {
  return goursatsurface(surface,p0,r,roots);
}
  
vec3 getnormal(Surface surface, vec3 p) {
  return goursatnormal(surface,p);
}
  
bool solve(Surface surface, vec3 p, vec3 r, float tmin, inout Result result) {
  vec4 roots;
  int nroots = dosurface(surface,p,r,roots);
  // Find smallest root greater than tmin.
  float t = result.t;
  for (int i = 0; i < 4; i++) {
    if (i == nroots) break;
    if (roots[i] > tmin && roots[i] < t) {
      t = roots[i];
    }
  }
  if (t == result.t) return false;
  p += t*r;
  vec3 n = getnormal(surface, p);
  if (dot(n,r) > 0.0) n = -n;
  n = normalize(n);
  vec3 basecolor = abs(n);
  if (surface.colorscheme == 1) {
    basecolor *= gridline(p);
  }
  result = Result(p,n,basecolor,t);
  return true;
}

// Interesting parameters from:
// https://www.mathcurve.com/surfaces.gb/goursat/goursat.shtml
vec4 goursatparams(int i) {
  int index = 0;
  if (i == index++) return vec4(0,-1,0,1);
  if (i == index++) return vec4(0,-2,2,1);
  if (i == index++) return vec4(-1,1,1,1);
  if (i == index++) return vec4(-1,-0.25,0.25,1);
  if (i == index++) return vec4(-0.5,-1,0.5,1);
  if (i == index++) return vec4(-0.5,1,-1.5,1);
  if (i == index++) return vec4(-1,4,-6,1);
  if (i == index++) return vec4(-1,1,1,1);
  if (i == index++) return vec4(-1,2,-2,1);
  else return vec4(-0.333,-0.666,0.666,1);
}

int nparams = 10;

int imod(int n, int m) {
    return n-n/m*m;
}

bool scene(vec3 p0, vec3 r, out vec3 col) {
  // Solve from closest point to origin.
  // This make p0.r = 0.
  float tmin = -dot(p0,r);
  p0 += tmin*r;
  Result res = Result(vec3(0),vec3(0),vec3(0),1e8);
  float ttime = 0.2*iTime;
  float rtime0 = floor(ttime);
  float rtime1 = floor(ttime+0.5);
  ttime = fract(2.0*ttime);
  vec4 params = vec4(0);
  int selection = int(iMouse.w)%(nparams+1);
  if (selection == 0) {
    params = mix(goursatparams(int(rtime0)%nparams),
                 goursatparams(int(rtime1)%nparams),
                 ttime);
  } else {
    params = goursatparams(selection-1);
  }
  Surface surface = Surface(params,vec3(0),1);
  if (!solve(surface,p0,r,-tmin,res)) return false;
  col = applylighting(res.basecolor,res.p,res.n,r);
  return true;
}

vec3 transform(vec3 p) {
  if (dorotate) {
    float t = iTime;
    p.yz = rotate(p.yz, 0.1*t);
    p.zx = rotate(p.zx, 0.222*t);
  }
  return p;
}

vec4 mainFun(vec3 p, vec3 r, vec3 rcentre) {
  p = transform(p);
  r = normalize(transform(r));
  rcentre = normalize(transform(rcentre));
  
  // Screenspace ray derivatives
  vec3 drdx = dFdx(r);
  vec3 drdy = dFdy(r);
  vec4 color = vec4(0);
  float k = dot(r,rcentre);
  // Just antialias the central part of the image.
  float aa0 = float(k > 0.96 ? 3 : k > 0.9 ? 2 : 1);
  float aa = aa0;
  if (iMouse.z > 0.0) aa = 1.0; // Just for comparison
  for (float i = 0.0; i < aa; i++) {
    for (float j = 0.0; j < aa; j++) {
      vec3 col1;
      // What to do with partially visible pixels?
      if (scene(p,normalize(r+(i-0.5*aa)/aa*drdx+(j-0.5*aa)/aa*drdy),col1)) {
        color += vec4(col1,1.0);
      }
    }
  }
  color /= aa*aa;
  if (dot(r,rcentre) > 0.999) color.b = 1.0;
  //color *= sqrt(aa0/3.0); // Show aa bands
  if (alert) color.x = 1.0;
  return pow(color,vec4(0.4545));
}

void mainVR(out vec4 fragColor, vec2 fragCoord, vec3 viewer, vec3 ray) {
  vec3 sceneoffset = vec3(0,0,8);
  // rcentre should be ray through centre of clipspace.
  vec3 rcentre = (inverse(iView)*vec4(0,0,-1,0)).xyz;
  fragColor = mainFun(viewer+sceneoffset,ray,rcentre);
}

vec3 mousetransform(vec3 p) {
  if (iMouse.x > 0.0) {
    float theta = (2.0*iMouse.y-iResolution.y)/iResolution.y*PI;
    float phi = (2.0*iMouse.x-iResolution.x)/iResolution.x*PI;
    p.yz = rotate(p.yz,theta);
    p.zx = rotate(p.zx,-phi);
  }
  return p;
}

void mainImage( out vec4 fragColor, vec2 fragCoord) {
  float scale = 1.0;
  float camera = 4.0;
  vec2 uv = scale*(2.0*fragCoord.xy - iResolution.xy)/iResolution.y;
  vec3 p = vec3(0,0,camera);
  vec3 r = vec3(uv,-2);
  vec3 rcentre = vec3(0,0,-1);
  p = mousetransform(p);
  r = mousetransform(r);
  rcentre = mousetransform(rcentre);
  fragColor = mainFun(p,r,rcentre);
}
