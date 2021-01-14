////////////////////////////////////////////////////////////////////////////////
//
// Created by Matthew Arcus, 2018
//
// A tetrahedron contains a parallelepiped with edges between the
// midpoints of its edges, and can also be inscribed uniquely in a
// parallepiped sharing its vertices.

// Mouse to move, 'e' and 'i' toggle display of internal and external
// parallelepipeds.
//
// Raytraced spheres and cylinders, solved analytically
//
////////////////////////////////////////////////////////////////////////////////


precision highp float;

bool dorotate = true;

const float PI = 3.14159;

struct Ray {
  vec3 q;               // origin
  vec3 d;               // direction
};

struct Hit {
  float t;      // solution to p=q+t*d
  vec3 n;       // (unnormalized) normal
  int id;       // what was hit
};

struct Sphere {
  float r2;      // radius squared
  vec3 p;       // centre
  int id;
};

struct Cylinder {
  float r2;      // radius squared
  // points s and t are end points of cylinder
  vec3 s,t;
  int id;
};

// Use vec2 return?
// Solve Ax^2 + 2Bx + C = 0
bool quadratic0(float A, float B, float C, out float x0, out float x1) {
   float D = B*B - A*C;
   if (D < 0.0) return false;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x0 = (-B-D)/A;
   x1 = C/(A*x0);
   if (B < 0.0) {
     // return smallest root first
     float t = x0; x0 = x1; x1 = t;
   }
   return true;
}

bool quadratic(float A, float B, float C, out float x) {
  float x1,x2;
  if (!quadratic0(A,B,C,x1,x2)) return false;
  if (x1 > 0.0) x = x1;
  else if (x2 > 0.0) x = x2;
  else return false;
  return true;
}

bool intersectSphere(Sphere s, Ray ray, inout Hit hit) {
  vec3 p = s.p;
  float r2 = s.r2;
  float c2 = dot(p,p);
  vec3 q = ray.q, d = ray.d;
  // |q + t*d - p|^2 = r^2
  float A = 1.0;
  float B = dot(q-p,d);
  float C = dot(q,q)-2.0*dot(q,p)+c2-r2;
  float t;
  if (!quadratic(A,B,C,t)) return false;
  if (t < 0.0 || t >= hit.t) return false;
  // Normal is the radial vector of sphere
  // We normalize it later
  hit = Hit(t, q+t*d-p, s.id);
  return true;
}

bool cylinder(vec3 q, vec3 s, float r2, vec3 p, vec3 r, out float t, out vec3 normal) {
  vec3 n = s-q; // Line direction
  float k2 = dot(n,n);
  vec3 p1 = p-q; // Move p to line space
  float rs = dot(r,n);
  float ps = dot(p1,n);
  float pr = dot(p1,r);
  float pp = dot(p1,p1);
  float A = 1.0 - rs*rs/k2;
  float B = pr - ps*rs/k2;
  float C = pp - ps*ps/k2 - r2;
  if (!quadratic(A,B,C,t)) {
    return false;
  } else {
    p1 += t*r; // Final point in line space
    float lambda = dot(p1,n)/k2;
    //if (lambda < 0.0 || lambda > 1.0) return false;
    normal = p1-lambda*n;
    return true;
  }
}

bool intersectCylinder(Cylinder c, Ray ray, inout Hit hit) {
  vec3 s0 = c.s;
  vec3 s1 = c.t;
  float r2 = c.r2;
  vec3 q = ray.q, d = ray.d;
  float t;
  vec3 normal;
  if (!cylinder(s0,s1,r2,q,d,t,normal) || t >= hit.t) {
    return false;
  } else {
    // Normal is the radial vector of cylinder
    hit = Hit(t, normal, c.id);
    return true;
  }
}

vec2 rotate(vec2 p, float t) {
  return p * cos(-t) + vec2(p.y, -p.x) * sin(-t);
}

vec3 transform(in vec3 p) {
  if (dorotate) {
    float t = iTime;
    p.zx = rotate(p.zx,t * 0.2);
  }
  if (iMouse.x > 0.0) {
    float theta = -(2.0*iMouse.y-iResolution.y)/iResolution.y*PI;
    float phi = -(2.0*iMouse.x-iResolution.x)/iResolution.x*PI;
    p.zx = rotate(p.zx,phi);
    p.yz = rotate(p.yz,-theta);
  }
  return p;
}

// (p+P)(q+Q) = pq + pQ + qP + PQ
vec4 qmul(vec4 p, vec4 q) {
  vec3 P = p.xyz, Q = q.xyz;
  return vec4(p.w*Q+q.w*P+cross(P,Q),p.w*q.w-dot(P,Q));
}

// Draw centroid at origin
const vec3 A = vec3(1,1,1);
const vec3 B = vec3(1,-1,-1);
const vec3 C = vec3(-1,1,-1);
const vec3 D = vec3(-1,-1,1);
const vec3 points[] = vec3[](A,B,C,D);

bool intersectScene(Ray ray, out Hit hit) {
  float t = 0.2*iTime;
  float s2 = 0.1*0.1;
  float c2 = 0.04*0.04;
  hit = Hit(1e8,vec3(0),0);
  for (int i = 0; i < 4; i++) {
    vec3 X = points[i];
    //intersectSphere(Sphere(s2,X,0),ray,hit);
    for (int j = i+1; j < 4; j++) {
      float a = 0.5+mod(0.5*iTime,3.0);
      vec3 Y = points[j];
      //intersectCylinder(Cylinder(c2,X,Y,1),ray,hit);
      intersectSphere(Sphere(s2,a*X+(1.0-a)*Y,3),ray,hit);
      intersectSphere(Sphere(s2,a*Y+(1.0-a)*X,3),ray,hit);
      for (int k = 0; k < 4; k++) {
        if (k == i || k == j) continue;
        vec3 Z = points[k];
        intersectCylinder(Cylinder(c2,a*X+(1.0-a)*Z,a*Y+(1.0-a)*Z,2),ray,hit);
        intersectCylinder(Cylinder(c2,a*Z+(1.0-a)*X,a*Z+(1.0-a)*Y,2),ray,hit);
      }
#if 0
      // Find points on X,Y distance k2 from origin
      float a = dot(X-Y,X-Y);
      float b = 2.0*dot(X-Y,Y);
      float c = dot(Y,Y) - k2;
      float d = b*b -4.0*a*c;
      if (d >= 0.0) {
        d = sqrt(d);
        float x0 = (-b-d)/(2.0*a);
        float x1 = (-b+d)/(2.0*a);
        intersectSphere(Sphere(s2,x0*X+(1.0-x0)*Y,1),ray,hit);
        intersectSphere(Sphere(s2,x1*X+(1.0-x1)*Y,1),ray,hit);
      }
#endif
    }
  }
  return hit.t < 1e8;
}

vec3 light = vec3(0);
float ambient = 0.2;
float diffuse = 0.8;

vec3 getColor(int i) {
  if (i == 0) return vec3(1,0,0);
  if (i == 1) return vec3(1,1,0);
  if (i == 2) return vec3(0,1,0);
  if (i == 3) return vec3(0,0,1);
  if (i == 4) return vec3(0,1,1);
  return vec3(1,1,1);
}

vec4 solve(Ray r) {
  Hit hit;
  if (!intersectScene(r,hit)) {
    return vec4(0);
  } else {
    vec3 n = normalize(hit.n);
    if (dot(r.d,n) > 0.0) n *= -1.0;
    vec3 baseColor = 0.7*getColor(hit.id);
    vec3 color = baseColor.xyz*(ambient+diffuse*max(0.0,dot(light,n)));
    float specular = pow(max(0.0,dot(reflect(light,n),r.d)),2.0);
    color += 0.5*specular*vec3(1.0,1.0,1.0);
    color *= clamp(1.0 - (hit.t-3.0)/5.0,0.0,1.0);
    return vec4(sqrt(color),1.0);
  }
}

void mainVR( out vec4 fragColor, in vec2 fragCoord, vec3 fragRayOrigin, vec3 fragRayDir) {
  float t = 0.5*iTime;
  light = normalize(vec3(0.5,1.0,-1.0));
  ambient = 0.5;
  diffuse = 1.0-ambient;
  fragColor = solve(Ray(fragRayOrigin+vec3(0,0,6),fragRayDir));
}


void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 uv = 2.0*fragCoord.xy/iResolution.xy - 1.0;
  vec3 p = vec3(0);
  // "screen" coordinate
  vec3 d = vec3(iResolution.x/iResolution.y * uv.x, uv.y, 2);
  p = transform(p);
  d = transform(d);
  d = normalize(d);
  mainVR(fragColor,fragCoord,p,d);
}
